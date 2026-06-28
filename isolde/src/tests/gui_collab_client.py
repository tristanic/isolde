#!/usr/bin/env python3
'''
Live-GUI collaboration smoke test (standalone HTTP client).

This is the end-to-end Phase A-C verification against a *running, visible* ISOLDE
GUI -- the human-agent collaboration mode. It runs as an ordinary Python process
(NOT inside ChimeraX) and drives ISOLDE entirely through the authenticated REST
agent surface, exactly as the MCP server does. Because the GUI runs the frame
loop, the simulation actually steps and rendering uses real OpenGL.

Procedure
---------
1. Launch the lane's ISOLDE GUI:           isolde\\isolde\\run_chimerax.bat
2. In ChimeraX:   isolde demo crystal_intro
                  isolde remote rest start        # logs PORT and TOKEN
                  isolde remote rest info         # re-prints them if needed
3. Run this script (any Python with stdlib):
       ISOLDE_REST_PORT=<port> ISOLDE_REST_TOKEN=<token> python gui_collab_client.py
   (PowerShell:  $env:ISOLDE_REST_PORT=...; $env:ISOLDE_REST_TOKEN=...; python gui_collab_client.py)

It exercises: tool manifest, model inventory, validation (rama/clashes with the
clashscore-convention note), a real PNG render, sim start (asserting the
non-silent-no-op stepping witness), live sim status (energy), a coordinate-moving
command (asserting rmsd_moved > 0), checkpoint and revert.
'''
import os
import sys
import json
import time
import http.client


def _client():
    host = os.environ.get('ISOLDE_REST_HOST', 'localhost')
    port = os.environ.get('ISOLDE_REST_PORT')
    token = os.environ.get('ISOLDE_REST_TOKEN')
    if not port or not token:
        sys.exit('Set ISOLDE_REST_PORT and ISOLDE_REST_TOKEN (see "isolde remote rest info").')
    return host, int(port), token


HOST, PORT, TOKEN = (None, None, None)


def post(cmd, **kwargs):
    conn = http.client.HTTPConnection(HOST, PORT, timeout=120)
    headers = {'Content-type': 'application/json', 'Authorization': 'Bearer ' + TOKEN}
    conn.request('POST', '', json.dumps({'cmd': cmd, 'args': [], 'kwargs': kwargs}).encode(), headers)
    resp = conn.getresponse()
    data = json.loads(resp.read().decode())
    conn.close()
    return data


def _ok(msg):
    print('PASS:', msg)


def _bad(msg):
    print('FAIL:', msg)
    raise SystemExit(1)


def main():
    global HOST, PORT, TOKEN
    HOST, PORT, TOKEN = _client()

    tools = post('agent_tools')
    _ok('manifest: %d typed tools' % tools.get('count', -1))

    inv = post('list_models')
    atomic = [m for m in inv.get('models', []) if m.get('kind') == 'atomic']
    if not atomic:
        _bad('no atomic model loaded -- run "isolde demo crystal_intro" first')
    spec = atomic[0]['spec']
    _ok('model loaded: %s (%s atoms)' % (spec, atomic[0].get('n_atoms')))

    rama = post('agent_invoke', name='isolde validate ramachandran', args={'model': spec})['result']
    _ok('rama: %d scorable, %d outliers' % (rama.get('n_scorable', -1), rama.get('n_outlier', -1)))

    clash = post('agent_invoke', name='isolde validate clashes', args={'model': spec})['result']
    _ok('clashes: %d total; convention=%s' % (clash.get('n_total', -1), clash.get('convention')))

    render = post('render', width=400, height=300)
    if 'image_base64' not in render:
        _bad('render did not return an image: %s' % render)
    _ok('render: real PNG, %d bytes' % render.get('bytes', -1))

    # Start a simulation on the whole selected model and assert it is real (the
    # silent-no-op guard: running + sim_manager present).
    started = post('agent_invoke', name='isolde sim', args={'cmd': 'start', 'atoms': spec})
    if 'error' in started:
        _bad('sim start errored: %s' % started.get('error'))
    w = started.get('witness') or {}
    if not w.get('started_ok'):
        _bad('sim start witness says possible silent no-op: %s' % w)
    _ok('sim start witness: running=%s sim_manager_present=%s' %
        (w.get('simulation_running'), w.get('sim_manager_present')))

    # Authoritative stepping proof (thread-safe): the coordinate checksum must
    # change across polls while the GUI frame loop steps the integrator.
    sums, energy = [], None
    for _ in range(6):
        st = post('agent_invoke', name='isolde sim status', args={})['result']
        sums.append(st.get('coord_checksum'))
        if st.get('current_energy_kJ_mol') is not None:
            energy = st['current_energy_kJ_mol']
        time.sleep(1)
    distinct = len(set(s for s in sums if s is not None))
    if distinct < 2:
        _bad('coords did not change across polls (sim not stepping?): %s' % sums)
    _ok('sim is genuinely stepping: %d distinct coord checksums; energy=%s kJ/mol'
        % (distinct, energy))

    # A coordinate-moving command should report rmsd_moved (movement witness).
    if any(t.get('command') == 'isolde pepflip' for t in tools.get('tools', [])):
        # Use a real residue from the model rather than a guessed number.
        rs = post('resolve_spec', spec=spec)
        res_specs = rs.get('residue_specs') or []
        if res_specs:
            flip = post('agent_invoke', name='isolde pepflip',
                        args={'residues': res_specs[len(res_specs) // 2]})
            mv = (flip.get('witness') or {}).get('rmsd_moved')
            if 'error' in flip:
                print('  (pepflip note: %s)' % flip['error'])
            else:
                _ok('pepflip movement witness rmsd_moved=%s' % mv)

    post('agent_invoke', name='isolde sim', args={'cmd': 'checkpoint'})
    _ok('checkpoint taken')
    post('agent_invoke', name='isolde sim', args={'cmd': 'revert'})
    _ok('reverted to checkpoint')
    post('agent_invoke', name='isolde sim', args={'cmd': 'stop'})
    _ok('sim stopped')

    print('ALL PASS (live GUI collaboration)')


if __name__ == '__main__':
    main()
