# ISOLDE MCP server

A standalone [Model Context Protocol](https://modelcontextprotocol.io) server that
exposes a running ISOLDE session as typed MCP tools, so an LLM agent (e.g. Claude)
can drive ISOLDE **alongside a human in a live GUI session** — the GUI is the
safety net, and the agent's actions are echoed into the ChimeraX log so the human
can follow (and, if need be, veto) every move (see *Human-in-the-loop*, below).

It is a **separate process** from ChimeraX and depends only on the `mcp` package
plus the standard library — it talks to ISOLDE exclusively over the authenticated
localhost REST API.

## 1. Start ISOLDE's REST server (inside ChimeraX)

```
isolde remote rest start
isolde remote rest info      # prints the port and bearer token
```

For long-running / unattended automation you may also pass `allowRun true` to
enable the (arbitrary-command) `run` method — off by default, as it is the one
acknowledged arbitrary-code path.

## 2. Run the MCP server (separate terminal)

```
pip install mcp
export ISOLDE_REST_HOST=localhost      # default
export ISOLDE_REST_PORT=<port>         # from 'isolde remote rest info'
export ISOLDE_REST_TOKEN=<token>       # from 'isolde remote rest info'
python path/to/chimerax/isolde/mcp/server.py
```

(On Windows PowerShell use `$env:ISOLDE_REST_PORT="..."` etc.)

Point your MCP client (Claude Desktop / Claude Code) at this command via its MCP
server configuration (stdio transport).

## What it exposes

**1. One tool per agent-safe ISOLDE command** (currently ~40), auto-discovered at
start-up from the REST manifest (the `agent_tools` method). Tool names are the
underscore form of the command — e.g. `isolde_sim`, `isolde_validate_ramachandran`,
`isolde_restrain_torsions`, `isolde_pepflip`, `isolde_stepto` — and their input
schemas are generated from each command's `CmdDesc` typing. **Exposing a command
in ISOLDE makes it appear here with no change to this server** (capture is
automatic; exposure is the opt-in `agent_safe` flag in `cmd/agent_manifest.py`).

**2. Perception / query helpers** (read-only structured state an eyeless agent
needs): `isolde_list_models`, `isolde_describe_model` (chains + one-letter
sequences + breaks + ligands), `isolde_resolve_spec` (confirm what an atom-spec
selects), `isolde_problem_zones` (clustered problem sites — "where do I focus?"),
`isolde_restraint_summary`, `isolde_residue_info` (rotamer / Rama / SS / B),
`isolde_map_info` (Rwork/Rfree + MDFF state — the density-fit query),
`isolde_bfactor_outliers`, `isolde_view_state`.

**3. Action helpers** that are not plain commands: `isolde_tug` (a transient
position-restraint spring — ISOLDE's gentle, geometry-preserving "move" primitive;
needs a running sim) and `isolde_set_mdff` (enable/disable or, rarely, re-weight
map coupling — with the anti-force-fit warning baked into its description).

**4. Async jobs** for long-running ops: `isolde_job_start` / `isolde_job_poll`
(used for `isolde brefine` / `brsr` B-factor refinement).

### Witnesses (silent-no-op guard)

Mutating commands return a **witness** alongside the result — a verifiable
statement of what actually changed (`rmsd_moved`, a restraint count delta, a sim
stepping/`stop_requested` signal) — so the agent can tell real work from a
command that succeeded without doing anything. Many ISOLDE effects are realized
*over subsequent simulation frames* (the integrator settles toward a target), so
witnesses are honest about deferred effects and tell the agent to confirm by
polling `isolde sim status` / `isolde_map_info`.

### Human-in-the-loop transparency

Because the agent shares a live GUI with a human, every **mutating** agent action
is echoed to the ChimeraX log as a concise `[agent] …` line plus its witnessed
outcome (e.g. `[agent] isolde restrain torsions #1.2  → +815 restraints`). Read-only
queries stay silent to avoid spamming the log. The full per-call log is always
returned to the agent regardless.

### Design philosophy (read this first)

The most important thing the server gives an agent is not a tool — it is ISOLDE's
**design philosophy**: ten principles describing what the target *is* and how to
use it without getting the user into a mess (the map is the evidence not the
scores; "non-outlier" ≠ "correct"; never force-fit; small errors are the dangerous
ones; let the physics settle; you work alongside a human; …). The single source of
truth is [`src/cmd/agent_philosophy.py`](../cmd/agent_philosophy.py), surfaced as:

- the MCP **`initialize` instructions** (the condensed constitution, always in the
  agent's context);
- the **`isolde://philosophy`** MCP resource (the full text, on demand);
- the **manifest preamble** (top of the `agent_tools` response); and
- **soft guardrails** — advisory `philosophy_flags` attached to a result when a
  move looks like it violates a principle (e.g. raising MDFF coupling, a very
  large single-step shift, mistaking a clean metric for correctness). Guardrails
  *flag and announce*, they never block — the GUI + human remain the safety net.

The standalone server fetches this over REST (`agent_philosophy`) at start-up, so
it stays in lockstep with ISOLDE; against an older ISOLDE without that method it
degrades gracefully (no instructions/resource, tools still work).

### Fit metrics, briefly

There is no single "make it fit" knob. For crystallographic data, `isolde_map_info`
reports Rwork/Rfree (which can rise slightly after a genuine improvement, pending
B-factor re-refinement). For cryo-EM there is no reliable scalar metric — visual
confirmation is primary. Raising MDFF coupling to chase fit force-fits atoms and
wrecks geometry; the aim is a low-energy model that *also* fits (principle 4 —
`set_mdff` will flag a coupling increase).

## Render (agent vision) — currently parked

An `isolde_render` tool (capturing the human's view plus orthogonal views as image
content) is **fully implemented but disabled**, pending a Qt 6.11 regression in
which offscreen image capture breaks the main window's resize-reshaping
(reproducible in vanilla ChimeraX 1.13.dev, absent in 1.12). While parked, the REST
`render` endpoint is not registered and the MCP tool is omitted from the manifest.
To re-enable once Qt is fixed: set `ENABLE_RENDER = True` in
`remote_control/rest_server/server.py`, uncomment the `isolde_render` entry in
`server.py`'s `_EXTRA_TOOLS`, and verify a capture no longer breaks resizing.

## Running environment

The MCP server works against any running ISOLDE that has the REST server started.
In practice that is a **live GUI session**: simulations and rendering are
frame-/GL-driven, so meaningful model-building work needs the GUI event loop
running (fully headless `--offscreen` operation is not available on all platforms,
including the primary Windows dev host, which lacks an offscreen GL backend).

## Status

Built and exercised end-to-end: the REST substrate (auth, auto-discovered typed
manifest, invoke, witnesses, perception helpers, async jobs) is implemented and
tested, and the surface has been driven over a real stdio MCP round-trip against a
live-GUI ISOLDE (an agent listing tools and running a perceive → act → verify
loop). Render is the one piece held back, behind the Qt 6.11 fix described above.
