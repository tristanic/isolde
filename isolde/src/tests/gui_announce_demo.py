#!/usr/bin/env python3
'''
Live demo of the GUI-log transparency policy (run against a live GUI ISOLDE).

Watch the ChimeraX *log* while this runs: every MUTATING agent action should
print a "[agent] ..." line (action) followed by a "  -> ..." outcome line, while
read-only QUERIES (validate / sim status / list_models) print NOTHING -- that is
the whole point of the policy.

Prereq -- launch the lane GUI with the bootstrap (loads the demo + starts the
REST server + writes the info file). From the bundle dir:

    set ISOLDE_REST_INFO_FILE=%TEMP%\isolde_rest_info.json
    run_chimerax.bat --script src/tests/gui_bootstrap.py

Then run this client (it reads the same info file):

    ISOLDE_REST_INFO_FILE=$TEMP/isolde_rest_info.json python src/tests/gui_announce_demo.py

This client makes no assertions; it narrates what it sends and what you should
see appear in the GUI log, pausing between steps so you can follow along.
'''
import os
import sys
import json
import time
import http.client

INFO_FILE = os.environ.get('ISOLDE_REST_INFO_FILE', 'rest_info.json')
HOST = os.environ.get('ISOLDE_REST_HOST', 'localhost')
PAUSE = float(os.environ.get('ISOLDE_DEMO_PAUSE', '2.0'))


def _load_info():
    for _ in range(120):
        if os.path.exists(INFO_FILE):
            try:
                with open(INFO_FILE) as f:
                    info = json.load(f)
                if info.get('port') and info.get('token'):
                    return info['port'], info['token']
            except Exception:
                pass
        time.sleep(0.5)
    sys.exit('No %s with port+token after 60s -- did the GUI bootstrap run?' % INFO_FILE)


PORT, TOKEN = _load_info()


def post(cmd, **kw):
    c = http.client.HTTPConnection(HOST, PORT, timeout=120)
    c.request('POST', '', json.dumps({'cmd': cmd, 'args': [], 'kwargs': kw}).encode(),
              {'Content-type': 'application/json', 'Authorization': 'Bearer ' + TOKEN})
    r = json.loads(c.getresponse().read().decode())
    c.close()
    return r


def step(note, expect_log, fn):
    '''Narrate one step: what we send, and what should appear in the GUI log.'''
    print('\n>>> %s' % note)
    if expect_log:
        print('    GUI log should show:  %s' % expect_log)
    else:
        print('    GUI log should show:  (nothing -- read-only)')
    result = fn()
    # Echo a compact view of the witness/result so the console mirrors the log.
    w = result.get('witness') if isinstance(result, dict) else None
    if w:
        print('    witness: %s' % {k: w[k] for k in (
            'rmsd_moved', 'max_atom_shift', 'restraint_total_delta',
            'restraint_enabled_delta', 'started_ok', 'stopped', 'launched')
            if k in w})
    elif isinstance(result, dict) and 'error' in result:
        print('    (error: %s)' % result['error'])
    time.sleep(PAUSE)
    return result


def wait_for_sim(running=True, timeout=20):
    deadline = time.time() + timeout
    while time.time() < deadline:
        r = post('agent_invoke', name='isolde sim status', args={})
        res = r.get('result') if isinstance(r, dict) else None
        if isinstance(res, dict) and bool(res.get('simulation_running')) == running:
            return True
        time.sleep(0.5)
    return False


def main():
    print('Connected to ISOLDE REST on port %d. Pause=%.1fs between steps.' % (PORT, PAUSE))
    print('=> Keep the ChimeraX LOG visible to watch the [agent] lines appear.')

    # Identify the loaded demo model.
    inv = post('list_models')
    atomic = [m for m in inv.get('models', []) if m.get('kind') == 'atomic']
    if not atomic:
        sys.exit('No atomic model loaded -- run "isolde demo crystal_intro" in the GUI first.')
    spec = next((m['spec'] for m in atomic if m.get('is_selected')), atomic[0]['spec'])
    print('Driving model %s' % spec)

    # 1. QUERY: validation -- should be SILENT in the GUI log.
    step('QUERY: validate ramachandran (read-only)',
         None,
         lambda: post('agent_invoke', name='isolde validate ramachandran', args={'model': spec}))

    # 2. MUTATING: start a simulation.
    step('MUTATING: start simulation',
         '[agent] isolde sim start %s   ->  simulation started' % spec,
         lambda: post('agent_invoke', name='isolde sim', args={'cmd': 'start', 'atoms': spec}))
    if not wait_for_sim(True):
        print('    (sim did not report running; tug/pepflip steps may no-op)')

    # 3. QUERY: sim status -- SILENT.
    step('QUERY: sim status (read-only)',
         None,
         lambda: post('agent_invoke', name='isolde sim status', args={}))

    # 4. MUTATING: impose torsion restraints on a chunk -> "+N restraints".
    step('MUTATING: restrain backbone torsions on a few residues',
         '[agent] isolde restrain torsions ...   ->  +N restraints',
         lambda: post('agent_invoke', name='isolde restrain torsions', args={'residues': spec}))

    # 5. MUTATING: tug a CA atom (helper that bypasses invoke_command).
    def do_tug():
        ri = post('residue_info', spec=spec)
        target = None
        # Tug the first residue's CA toward a point 2 A away, via to_spec-free target.
        rs = post('resolve_spec', spec=spec + '@CA')
        if rs.get('valid') and rs.get('centroid'):
            c = rs['centroid']
            target = [c[0] + 2.0, c[1], c[2]]
        return post('tug', spec=spec + '@CA', target=target)
    step('MUTATING: tug all CA atoms ~2 A (transient spring)',
         '[agent] tug %s@CA   ->  ~2.0 A (N atoms)' % spec,
         do_tug)

    # 6. MUTATING: a peptide flip on one residue.
    flip_spec = spec + '/A:78' if '/' not in spec else spec
    step('MUTATING: pepflip one residue',
         '[agent] isolde pepflip %s   ->  moved / target set' % flip_spec,
         lambda: post('agent_invoke', name='isolde pepflip', args={'atoms': flip_spec}))

    # 7. MUTATING: release the torsion restraints -> "released".
    step('MUTATING: release torsion restraints',
         '[agent] isolde release torsions ...   ->  N restraints released',
         lambda: post('agent_invoke', name='isolde release torsions', args={'residues': spec}))

    # 8. MUTATING: stop the simulation.
    step('MUTATING: stop simulation',
         '[agent] isolde sim stop   ->  simulation stopped',
         lambda: post('agent_invoke', name='isolde sim', args={'cmd': 'stop'}))

    print('\nDemo complete. The GUI log should contain a [agent] line for every '
          'MUTATING step above, and none for the two QUERY steps.')


if __name__ == '__main__':
    main()
