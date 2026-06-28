# @Author: Tristan Croll
# @Date:   27-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 27-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Standalone MCP server fronting a running ISOLDE session.

On start-up it fetches the auto-generated tool manifest from ISOLDE's REST API
(the ``agent_tools`` method), then emits one MCP tool per agent-safe ISOLDE
command plus the helper tools defined below: perception/query helpers (model
inventory, model description, atom-spec resolution, problem zones, restraint /
residue / map info, B-factor outliers, view state), action helpers (tug, MDFF
coupling), and async-job control. Tool calls are forwarded to the authenticated
REST endpoints. (An ``isolde_render`` image tool is implemented but PARKED behind
a Qt 6.11 regression -- see _EXTRA_TOOLS and README.md.)

Run it as its own process (NOT inside ChimeraX) -- it depends only on the
``mcp`` package and the standard library:

    pip install mcp
    export ISOLDE_REST_PORT=<port>      # from 'isolde remote rest info'
    export ISOLDE_REST_TOKEN=<token>    # from 'isolde remote rest info'
    python server.py

See README.md. The surface has been driven over a real stdio MCP round-trip
against a live-GUI ISOLDE; render is the one piece held back pending the Qt fix.
'''
import os
import sys
import json
import http.client


# --------------------------------------------------------------------------
# Minimal authenticated REST client (sync; the MCP layer is low-traffic).
# --------------------------------------------------------------------------
class IsoldeREST:
    def __init__(self, host, port, token, timeout=120):
        self.host, self.port, self.token, self.timeout = host, port, token, timeout

    def _call(self, method, body=None):
        conn = http.client.HTTPConnection(self.host, self.port, timeout=self.timeout)
        headers = {'Content-type': 'application/json',
                   'Authorization': 'Bearer ' + self.token}
        payload = json.dumps(body).encode('utf-8') if body is not None else None
        conn.request(method, '', payload, headers)
        resp = conn.getresponse()
        data = resp.read().decode('utf-8')
        conn.close()
        if resp.status == 401:
            raise PermissionError('ISOLDE REST rejected the bearer token (401)')
        try:
            return json.loads(data)
        except Exception:
            return {'error': 'non-JSON response', 'raw': data, 'status': resp.status}

    def post(self, cmd, **kwargs):
        return self._call('POST', {'cmd': cmd, 'args': [], 'kwargs': kwargs})

    def tools(self):
        return self.post('agent_tools').get('tools', [])

    def philosophy(self):
        '''Fetch ISOLDE's design philosophy (instructions + full principles). Returns
        {} against an older ISOLDE that lacks the agent_philosophy method, so the
        caller degrades gracefully.'''
        try:
            result = self.post('agent_philosophy')
        except Exception:
            return {}
        if not isinstance(result, dict) or 'error' in result or 'instructions' not in result:
            return {}
        return result

    def invoke(self, command, args):
        return self.post('agent_invoke', name=command, args=args)


# --------------------------------------------------------------------------
# Transport tools that are not ISOLDE commands (explicitly described here).
# --------------------------------------------------------------------------
_EXTRA_TOOLS = [
    {'name': 'isolde_list_models', 'rest': 'list_models',
     'description': 'List open models (atomic structures and maps) with their IDs, '
                    'names, chains and atom counts, so you can discover what is loaded.',
     'input_schema': {'type': 'object', 'properties': {}, 'additionalProperties': False}},
    {'name': 'isolde_describe_model', 'rest': 'describe_model',
     'description': 'Structured description of a model: chains with one-letter '
                    'sequences and modelled/unmodelled state, chain breaks, '
                    'ligands/non-standard residues, and water count. Use to '
                    'understand a model before building into it.',
     'input_schema': {'type': 'object', 'properties': {
         'model': {'type': 'string', 'description': "model id or spec, e.g. '#1.2'"}},
         'required': ['model'], 'additionalProperties': False}},
    {'name': 'isolde_problem_zones', 'rest': 'problem_zones',
     'description': 'Cluster a model\'s geometry/restraint problems (unsatisfied '
                    'restraints, rotamer/backbone/clash outliers) into spatial '
                    'zones, so you know where to focus. Returns zones with '
                    'centroid + member residues + a problem-type breakdown.',
     'input_schema': {'type': 'object', 'properties': {
         'model': {'type': 'string'},
         'outliers_only': {'type': 'boolean'}},
         'required': ['model'], 'additionalProperties': False}},
    {'name': 'isolde_restraint_summary', 'rest': 'restraint_summary',
     'description': 'Per restraint-type counts (total/enabled/unsatisfied) for a '
                    'model: what restraints are active and whether satisfied.',
     'input_schema': {'type': 'object', 'properties': {'model': {'type': 'string'}},
         'required': ['model'], 'additionalProperties': False}},
    {'name': 'isolde_residue_info', 'rest': 'residue_info',
     'description': 'Per-residue detail for targeted fixes: rotamer score/options, '
                    'Ramachandran case/score/phi-psi, secondary structure, B-factor '
                    'stats. Pass an atom-spec selecting one or more residues.',
     'input_schema': {'type': 'object', 'properties': {
         'spec': {'type': 'string', 'description': "e.g. '#1.2/A:58'"}},
         'required': ['spec'], 'additionalProperties': False}},
    {'name': 'isolde_map_info', 'rest': 'map_info',
     'description': 'Density/map-fit query: maps present, Rwork/Rfree '
                    '(crystallographic), and MDFF coupling state. For cryo-EM there '
                    'is no reliable scalar fit metric — use isolde_render for '
                    'visual confirmation. See the returned fit_guidance.',
     'input_schema': {'type': 'object', 'properties': {'model': {'type': 'string'}},
         'required': ['model'], 'additionalProperties': False}},
    {'name': 'isolde_set_mdff', 'rest': 'set_mdff',
     'description': 'Enable/disable or (rarely) re-weight MDFF map coupling. WARNING: '
                    'raising the coupling weight to force better fit ruins geometry; '
                    'ISOLDE\'s aim is a low-energy model that also fits. Prefer fixing '
                    'the model over cranking the weight.',
     'input_schema': {'type': 'object', 'properties': {
         'model': {'type': 'string'}, 'volume': {'type': 'string'},
         'enabled': {'type': 'boolean'},
         'coupling_constant': {'type': 'number'}},
         'required': ['model'], 'additionalProperties': False}},
    {'name': 'isolde_bfactor_outliers', 'rest': 'bfactor_outliers',
     'description': 'Atoms with extreme B-factors (a cheap problem sniffer, esp. '
                    'after refinement; caveat: genuinely mobile sidechains also '
                    'show high B).',
     'input_schema': {'type': 'object', 'properties': {
         'model': {'type': 'string'}, 'z_threshold': {'type': 'number'}},
         'required': ['model'], 'additionalProperties': False}},
    {'name': 'isolde_tug', 'rest': 'tug',
     'description': 'Tug atoms toward a target with a transient spring in a running '
                    'sim — ISOLDE\'s gentle, geometry-preserving way to move things: '
                    'express intent ("this belongs ~there") and let MD relax around '
                    'it, rather than hard-setting coordinates. Nudge a residue into '
                    'density or out of a clash, then watch it settle. Multi-atom '
                    'selections translate rigidly. Requires a running simulation; '
                    'release with release=true. Avoid strong spring_constant (force-fits).',
     'input_schema': {'type': 'object', 'properties': {
         'spec': {'type': 'string', 'description': 'atoms to tug'},
         'target': {'type': 'array', 'items': {'type': 'number'},
                    'description': '[x,y,z] in Angstroms to pull the centroid to'},
         'to_spec': {'type': 'string', 'description': 'alternatively, tug toward this spec\'s centroid'},
         'spring_constant': {'type': 'number'},
         'release': {'type': 'boolean'}},
         'required': ['spec'], 'additionalProperties': False}},
    {'name': 'isolde_view_state', 'rest': 'view_state',
     'description': 'Current selection + view orientation + ISOLDE stepper cursor.',
     'input_schema': {'type': 'object', 'properties': {},
         'additionalProperties': False}},
    {'name': 'isolde_resolve_spec', 'rest': 'resolve_spec',
     'description': 'Resolve a ChimeraX atom-spec and describe what it selects '
                    '(counts, bounding box, representative residues) before acting on it.',
     'input_schema': {'type': 'object',
                      'properties': {'spec': {'type': 'string',
                          'description': "ChimeraX atom-spec, e.g. '#1/A:58'"}},
                      'required': ['spec'], 'additionalProperties': False}},
    # isolde_render is PARKED pending a Qt 6.11 fix (offscreen capture breaks the
    # main window's resize-reshaping; see render.py / server.py ENABLE_RENDER).
    # The REST `render` endpoint is not registered while parked, so this tool is
    # omitted from the manifest. Re-enable together with server.py ENABLE_RENDER:
    # {'name': 'isolde_render', 'rest': 'render',
    #  'description': 'Render the human view + orthogonal views as PNG image content.',
    #  'input_schema': {'type': 'object', 'properties': {
    #      'width': {'type': 'integer', 'minimum': 1},
    #      'height': {'type': 'integer', 'minimum': 1},
    #      'extra_angles': {'type': 'array', 'items': {'type': 'number'}}},
    #      'additionalProperties': False}},
    {'name': 'isolde_job_start', 'rest': 'job_start',
     'description': 'Start a long-running ISOLDE command as a background job; '
                    'returns a job_id to poll with isolde_job_poll.',
     'input_schema': {'type': 'object', 'properties': {
         'name': {'type': 'string', 'description': 'command name, e.g. isolde brefine'},
         'args': {'type': 'object'}}, 'required': ['name'],
         'additionalProperties': False}},
    {'name': 'isolde_job_poll', 'rest': 'job_poll',
     'description': 'Poll the state (and result when finished) of a background job.',
     'input_schema': {'type': 'object',
                      'properties': {'job_id': {'type': 'string'}},
                      'required': ['job_id'], 'additionalProperties': False}},
]
_EXTRA_BY_NAME = {t['name']: t for t in _EXTRA_TOOLS}


_PHILOSOPHY_URI = 'isolde://philosophy'


def build_server(rest: 'IsoldeREST', philosophy=None):
    '''Construct the MCP server. Imported lazily so this file is importable
    (e.g. for inspection) even where the mcp package is absent.

    *philosophy* is the dict fetched from ISOLDE's agent_philosophy method (or {}
    against an older ISOLDE); when it carries a document, it is exposed as the
    ``isolde://philosophy`` resource.'''
    from mcp.server import Server
    from mcp.types import Tool, TextContent, ImageContent, Resource

    server = Server('isolde')
    philosophy = philosophy or {}
    philosophy_doc = philosophy.get('document')

    if philosophy_doc:
        @server.list_resources()
        async def list_resources():
            return [Resource(uri=_PHILOSOPHY_URI,
                             name='ISOLDE design philosophy',
                             description='Principles for driving ISOLDE safely '
                                         '(read before mutating a model).',
                             mimeType='text/plain')]

        @server.read_resource()
        async def read_resource(uri):
            if str(uri) == _PHILOSOPHY_URI:
                return philosophy_doc
            raise ValueError('unknown resource: %s' % uri)

    # Maps an MCP tool name to the canonical ISOLDE command name. Tool names are
    # MCP-sanitised (e.g. 'isolde ~ignore' -> 'isolde_not_ignore'), so they cannot
    # be un-munged back to the command by string surgery -- we carry the command
    # explicitly from the manifest. Populated by list_tools, read by call_tool.
    command_by_tool = {}

    @server.list_tools()
    async def list_tools():
        tools = []
        # One MCP tool per agent-safe ISOLDE command (auto-discovered manifest).
        for t in rest.tools():
            if 'error' in t:
                continue
            command_by_tool[t['name']] = t.get('command', t['name'].replace('_', ' '))
            tools.append(Tool(name=t['name'],
                              description=t.get('description', t['name']),
                              inputSchema=t.get('input_schema', {'type': 'object'})))
        # Plus the transport tools.
        for t in _EXTRA_TOOLS:
            tools.append(Tool(name=t['name'], description=t['description'],
                              inputSchema=t['input_schema']))
        return tools

    @server.call_tool()
    async def call_tool(name, arguments):
        arguments = arguments or {}
        # Transport tools forward to their REST endpoint directly.
        if name in _EXTRA_BY_NAME:
            rest_name = _EXTRA_BY_NAME[name]['rest']
            result = rest.post(rest_name, **arguments)
            if name == 'isolde_render' and isinstance(result, dict) \
                    and 'image_base64' in result:
                return [ImageContent(type='image', mimeType='image/png',
                                     data=result['image_base64'])]
            return [TextContent(type='text', text=json.dumps(result, indent=2))]
        # Otherwise it is an ISOLDE command tool: forward to agent_invoke using
        # the canonical command name from the manifest (the tool name is an
        # MCP-sanitised label and may not un-munge cleanly). Fall back to the
        # name if call_tool somehow precedes list_tools.
        command = command_by_tool.get(name) or name.replace('_', ' ')
        result = rest.invoke(command, arguments)
        return [TextContent(type='text', text=json.dumps(result, indent=2))]

    return server


def _config_from_env():
    host = os.environ.get('ISOLDE_REST_HOST', 'localhost')
    port = os.environ.get('ISOLDE_REST_PORT')
    token = os.environ.get('ISOLDE_REST_TOKEN')
    if not port or not token:
        sys.stderr.write(
            'ISOLDE MCP server: set ISOLDE_REST_PORT and ISOLDE_REST_TOKEN '
            '(see "isolde remote rest info" in ChimeraX).\n')
        raise SystemExit(2)
    return host, int(port), token


def main():
    import asyncio
    from mcp.server.stdio import stdio_server
    from mcp.server import InitializationOptions, NotificationOptions

    host, port, token = _config_from_env()
    rest = IsoldeREST(host, port, token)
    philosophy = rest.philosophy()      # {} against an older ISOLDE
    server = build_server(rest, philosophy)

    init_kwargs = dict(
        server_name='isolde',
        server_version='0.1.0',
        capabilities=server.get_capabilities(
            notification_options=NotificationOptions(),
            experimental_capabilities={}),
    )
    # Surface the philosophy as the initialize `instructions` (always in the
    # agent's context). Guard for older mcp SDKs whose InitializationOptions has
    # no such field -- the resource + manifest preamble still carry it.
    instructions = philosophy.get('instructions')
    if instructions and _supports_instructions(InitializationOptions):
        init_kwargs['instructions'] = instructions

    async def _run():
        async with stdio_server() as (read, write):
            await server.run(read, write, InitializationOptions(**init_kwargs))

    asyncio.run(_run())


def _supports_instructions(InitializationOptions):
    '''True if this mcp SDK's InitializationOptions accepts an `instructions` field.'''
    fields = getattr(InitializationOptions, 'model_fields', None) \
        or getattr(InitializationOptions, '__fields__', {})
    return 'instructions' in fields


if __name__ == '__main__':
    main()
