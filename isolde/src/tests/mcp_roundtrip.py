#!/usr/bin/env python3
'''
End-to-end MCP round-trip test (the deferred Phase D check).

Acts as a real MCP *client*: launches src/mcp/server.py over stdio, does the MCP
initialize handshake, lists tools, and calls a representative set — driving the
live ISOLDE GUI through MCP exactly as an LLM agent (Claude) would. Run with the
ChimeraX python (which has the mcp SDK):

    ISOLDE_REST_PORT=<port> ISOLDE_REST_TOKEN=<token> \
        "C:/Program Files/ChimeraX-Daily/bin/python.exe" mcp_roundtrip.py
'''
import os
import sys
import asyncio

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

SERVER = os.path.join(os.path.dirname(__file__), '..', 'mcp', 'server.py')


def _text(result):
    parts = []
    for c in result.content:
        if getattr(c, 'type', None) == 'text':
            parts.append(c.text)
        elif getattr(c, 'type', None) == 'image':
            parts.append('<image %s, %d b64 chars>' % (c.mimeType, len(c.data)))
    return '\n'.join(parts)


async def main():
    env = dict(os.environ)
    params = StdioServerParameters(command=sys.executable, args=[SERVER], env=env)
    print('>>> launching ISOLDE MCP server over stdio...')
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as session:
            info = await session.initialize()
            print('>>> MCP handshake OK — server:', info.serverInfo.name,
                  info.serverInfo.version)

            tools = (await session.list_tools()).tools
            print('>>> %d tools advertised. sample:' % len(tools))
            for t in tools[:12]:
                print('      -', t.name)

            async def call(name, **args):
                print('\n>>> call_tool %s %s' % (name, args))
                r = await session.call_tool(name, args)
                out = _text(r)
                print(out[:900])
                return out

            # --- the agentic perceive -> act -> verify loop, over MCP ---
            await call('isolde_list_models')
            await call('isolde_describe_model', model='#1.2')
            await call('isolde_validate_ramachandran', model='#1.2')        # command tool
            await call('isolde_map_info', model='#1.2')             # Rwork/Rfree
            await call('isolde_render', width=400, height=300)      # -> image content

            # act: start a sim, tug a residue, confirm, release, stop
            await call('isolde_sim', cmd='start', atoms='#1.2')
            await asyncio.sleep(2)
            spec = '#1.2/A:100'
            import json
            rs = await call('isolde_resolve_spec', spec=spec + '@CA')
            try:
                c = json.loads(rs).get('centroid')
            except Exception:
                c = None
            if c:
                await call('isolde_tug', spec=spec, target=[c[0] + 3.0, c[1], c[2]])
                await asyncio.sleep(4)
                await call('isolde_sim_status')                     # command tool
                await call('isolde_tug', spec=spec, release=True)
            await call('isolde_sim', cmd='stop')

            print('\n>>> MCP round-trip complete.')


if __name__ == '__main__':
    asyncio.run(main())
