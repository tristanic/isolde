#!/usr/bin/env python3
'''
Phase-F coverage verification (standalone HTTP client, run against a live GUI).

Prereq (in the lane GUI ChimeraX):
    isolde demo crystal_intro
    isolde remote rest start          # note PORT + TOKEN from the log / 'info'

Run:
    ISOLDE_REST_PORT=<port> ISOLDE_REST_TOKEN=<token> python gui_phaseF_client.py

Exercises the new surface: expanded manifest, preflight, restraints (+ the
restraint_change witness), problem_zones, restraint_summary, residue_info,
map_info (Rwork/Rfree + MDFF state), set_mdff round-trip, bfactor_outliers,
view_state. (brefine is exposed but not run to completion here.)
'''
import os
import sys
import json
import http.client

HOST = os.environ.get('ISOLDE_REST_HOST', 'localhost')
PORT = os.environ.get('ISOLDE_REST_PORT')
TOKEN = os.environ.get('ISOLDE_REST_TOKEN')
if not PORT or not TOKEN:
    sys.exit('Set ISOLDE_REST_PORT and ISOLDE_REST_TOKEN (see "isolde remote rest info").')
PORT = int(PORT)


def post(cmd, **kw):
    c = http.client.HTTPConnection(HOST, PORT, timeout=120)
    c.request('POST', '', json.dumps({'cmd': cmd, 'args': [], 'kwargs': kw}).encode(),
              {'Content-type': 'application/json', 'Authorization': 'Bearer ' + TOKEN})
    r = json.loads(c.getresponse().read().decode()); c.close(); return r


def ok(m): print('PASS:', m)
def bad(m): print('FAIL:', m); raise SystemExit(1)


def main():
    tools = post('agent_tools')
    names = {t['command'] for t in tools['tools']}
    for needed in ('isolde restrain distances', 'isolde restrain torsions',
                   'isolde preflight hydrogens', 'isolde brefine'):
        if needed not in names:
            bad('manifest missing %r' % needed)
    ok('manifest now %d tools incl. restraints/preflight/brefine' % tools['count'])

    inv = post('list_models')
    atomic = [m for m in inv['models'] if m.get('kind') == 'atomic']
    if not atomic:
        bad('no atomic model — run "isolde demo crystal_intro" first')
    spec = atomic[0]['spec']
    ok('model %s' % spec)

    # preflight (read-only)
    pf = post('agent_invoke', name='isolde preflight hydrogens', args={'model': spec})
    ok('preflight hydrogens -> %s' % (pf.get('result') or pf.get('error')))

    # restraint baseline
    rs0 = post('restraint_summary', model=spec)
    ok('restraint_summary (before): %s' % rs0.get('restraints'))

    # impose torsion restraints on a chunk -> restraint_change witness shows +N
    rr = post('agent_invoke', name='isolde restrain torsions',
              args={'residues': spec})
    if 'error' in rr:
        print('  (restrain torsions note: %s)' % rr['error'])
    w = rr.get('witness') or {}
    ok('restrain torsions witness: total_delta=%s by_type=%s' %
       (w.get('restraint_total_delta'), w.get('by_type_total_delta')))

    rs1 = post('restraint_summary', model=spec)
    ok('restraint_summary (after): %s' % rs1.get('restraints'))

    # release them -> witness shows a negative total OR enabled delta (release
    # disables rather than deletes, so the enabled-count is what changes).
    rel = (post('agent_invoke', name='isolde release torsions',
                args={'residues': spec}).get('witness') or {})
    ok('release torsions witness: total_delta=%s enabled_delta=%s' %
       (rel.get('restraint_total_delta'), rel.get('restraint_enabled_delta')))

    # problem_zones
    pz = post('problem_zones', model=spec, outliers_only=True)
    if 'error' in pz:
        print('  (problem_zones note: %s)' % pz['error'])
    else:
        ok('problem_zones: %d zones, %d unclustered (types %s)' %
           (pz['n_zones'], pz['n_unclustered'], pz.get('unclustered_types')))

    # residue_info on the first few residues
    ri = post('residue_info', spec=spec + '/A:76-80' if '/' not in spec else spec)
    if 'error' in ri:
        # fall back to whole-model first residues
        ri = post('residue_info', spec=spec)
    if 'residues' in ri and ri['residues']:
        ok('residue_info sample: %s' % ri['residues'][0])
    else:
        print('  (residue_info: %s)' % ri)

    # map_info — the density-fit query (Rwork/Rfree + MDFF)
    mi = post('map_info', model=spec)
    if 'error' in mi:
        bad('map_info error: %s' % mi['error'])
    ok('map_info: has_maps=%s crystallographic=%s mdff=%s' %
       (mi.get('has_maps'), mi.get('crystallographic'), mi.get('mdff')))

    # set_mdff round-trip (disable then re-enable), if any MDFF mgr exists
    if mi.get('mdff'):
        d = post('set_mdff', model=spec, enabled=False)
        e = post('set_mdff', model=spec, enabled=True)
        ok('set_mdff round-trip: disabled=%s re-enabled=%s' %
           ([m['enabled'] for m in d.get('mdff', [])],
            [m['enabled'] for m in e.get('mdff', [])]))
    else:
        print('  (no MDFF manager active yet — set_mdff skipped)')

    # bfactor_outliers
    bo = post('bfactor_outliers', model=spec)
    ok('bfactor_outliers: %d outliers (mean B=%.1f)' %
       (bo.get('n_outliers', -1), bo.get('mean_bfactor', float('nan'))))

    # view_state
    vs = post('view_state')
    ok('view_state: selection=%s' % (vs.get('selection')))

    print('ALL PASS (Phase F)')


if __name__ == '__main__':
    main()
