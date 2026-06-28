#!/usr/bin/env python3
'''
Control experiment: exercise the FULL agent interface against a live GUI with
NO render/capture calls at all (so neither the offscreen capture nor the
post-render resync nudge ever runs). If window-resize reshaping breaks after
this, the capture is exonerated and the trigger is elsewhere.

    ISOLDE_REST_PORT=<port> ISOLDE_REST_TOKEN=<token> python gui_norender_client.py
'''
import os, sys, json, http.client, time

PORT = os.environ['ISOLDE_REST_PORT']; TOKEN = os.environ['ISOLDE_REST_TOKEN']
PORT = int(PORT)


def post(cmd, **kw):
    c = http.client.HTTPConnection('localhost', PORT, timeout=60)
    c.request('POST', '', json.dumps({'cmd': cmd, 'args': [], 'kwargs': kw}).encode(),
              {'Content-type': 'application/json', 'Authorization': 'Bearer ' + TOKEN})
    r = json.loads(c.getresponse().read().decode()); c.close(); return r


def step(label, val):
    print('PASS: %s -> %s' % (label, (str(val)[:120])))


SPEC = '#1.2'


def main():
    step('agent_tools count', post('agent_tools').get('count'))
    inv = post('list_models')
    step('list_models', '%d models' % inv.get('count', -1))
    step('describe_model', 'chains=%s' % [c['chain_id'] for c in post('describe_model', model=SPEC).get('chains', [])])
    step('validate ramachandran', post('agent_invoke', name='isolde validate ramachandran', args={'model': SPEC})['result'].get('n_outlier'))
    step('validate clashes', post('agent_invoke', name='isolde validate clashes', args={'model': SPEC})['result'].get('n_total'))
    mi = post('map_info', model=SPEC)
    step('map_info Rwork/Rfree', [ (x.get('rwork'), x.get('rfree')) for x in mi.get('crystallographic', []) ])
    step('restraint_summary', post('restraint_summary', model=SPEC).get('restraints'))
    pz = post('problem_zones', model=SPEC)
    step('problem_zones', 'zones=%s unclustered=%s' % (pz.get('n_zones'), pz.get('n_unclustered')))
    step('residue_info', post('residue_info', spec=SPEC + '/A:100').get('residues', [{}])[0].get('spec'))
    step('preflight hydrogens', post('agent_invoke', name='isolde preflight hydrogens', args={'model': SPEC})['result'].get('status'))
    step('bfactor_outliers', post('bfactor_outliers', model=SPEC).get('n_outliers'))
    step('view_state', post('view_state').get('selection', {}).get('n_atoms'))

    # action: sim + tug + restraints + mdff  (still NO render)
    step('sim start', (post('agent_invoke', name='isolde sim', args={'cmd': 'start', 'atoms': SPEC}).get('witness') or {}).get('started_ok'))
    time.sleep(2)
    sums = [post('agent_invoke', name='isolde sim status', args={})['result'].get('coord_checksum') for _ in range(3)]
    step('sim stepping (distinct checksums)', len(set(sums)))
    c = post('resolve_spec', spec=SPEC + '/A:100@CA').get('centroid')
    if c:
        step('tug', post('tug', spec=SPEC + '/A:100', target=[c[0] + 3, c[1], c[2]]).get('n_tugged'))
        time.sleep(2)
        post('tug', spec=SPEC + '/A:100', release=True)
    step('restrain torsions', (post('agent_invoke', name='isolde restrain torsions', args={'residues': SPEC}).get('witness') or {}).get('restraint_total_delta'))
    post('agent_invoke', name='isolde release torsions', args={'residues': SPEC})
    if mi.get('mdff'):
        post('set_mdff', model=SPEC, enabled=False)
        step('set_mdff round-trip', [m['enabled'] for m in post('set_mdff', model=SPEC, enabled=True).get('mdff', [])])
    post('agent_invoke', name='isolde sim', args={'cmd': 'stop'})
    print('ALL DONE (no render performed)')


if __name__ == '__main__':
    main()
