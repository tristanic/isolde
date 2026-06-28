#!/usr/bin/env python3
'''
Heavy end-to-end workout of the agent interface against a live GUI, for a
post-run GUI-health check (resize/dock). Exercises perception, render (expected
to be a safe deferred no-op now), and multiple simulation / tug / restraint /
MDFF / refinement cycles.

    ISOLDE_REST_PORT=<port> ISOLDE_REST_TOKEN=<token> python gui_workout_client.py
'''
import os, sys, json, http.client, time

PORT = int(os.environ['ISOLDE_REST_PORT']); TOKEN = os.environ['ISOLDE_REST_TOKEN']
SPEC = '#1.2'
nfail = 0


def post(cmd, **kw):
    c = http.client.HTTPConnection('localhost', PORT, timeout=60)
    c.request('POST', '', json.dumps({'cmd': cmd, 'args': [], 'kwargs': kw}).encode(),
              {'Content-type': 'application/json', 'Authorization': 'Bearer ' + TOKEN})
    r = json.loads(c.getresponse().read().decode()); c.close(); return r


def inv(name, **args):
    return post('agent_invoke', name=name, args=args)


def show(label, val):
    print('  %-34s %s' % (label, str(val)[:110]))


def main():
    global nfail
    print('== discovery + perception ==')
    show('agent_tools', post('agent_tools').get('count'))
    show('list_models', post('list_models').get('count'))
    show('describe_model', [c['chain_id'] for c in post('describe_model', model=SPEC).get('chains', [])])
    for v in ('rama', 'rotamers', 'clashes', 'peptidebonds'):
        r = inv('isolde validate ' + v, model=SPEC).get('result', {})
        show('validate ' + v, {k: r[k] for k in r if k.startswith('n_')})
    for p in ('hydrogens', 'parameters', 'disulfides', 'altlocs'):
        show('preflight ' + p, inv('isolde preflight ' + p, model=SPEC).get('result', {}).get('status', '(ok)'))
    show('map_info', [(x.get('rwork'), x.get('rfree')) for x in post('map_info', model=SPEC).get('crystallographic', [])])
    show('problem_zones', 'zones=%s' % post('problem_zones', model=SPEC).get('n_zones'))
    show('residue_info', post('residue_info', spec=SPEC + '/A:84').get('residues', [{}])[0].get('rama'))
    show('bfactor_outliers', post('bfactor_outliers', model=SPEC).get('n_outliers'))
    show('view_state', post('view_state').get('selection', {}).get('n_atoms'))

    print('== render (expect safe deferred no-op) ==')
    rr = post('render', width=400, height=300)
    if rr.get('deferred'):
        show('render', 'deferred no-op OK')
    elif 'image_base64' in rr:
        show('render', 'returned image (%d bytes)?!' % rr.get('bytes', -1))
    else:
        show('render', rr.get('error'))

    print('== simulation / tug / restraint cycles ==')
    for cycle in range(2):
        print('  -- cycle %d --' % cycle)
        w = (inv('isolde sim', cmd='start', atoms=SPEC).get('witness') or {})
        show('sim start started_ok', w.get('started_ok'))
        if not w.get('started_ok'):
            nfail += 1
        time.sleep(2)
        sums = [inv('isolde sim status').get('result', {}).get('coord_checksum') for _ in range(3)]
        show('stepping (distinct checksums)', len(set(sums)))
        # tug a couple of residues
        for resnum in (60, 120):
            c = post('resolve_spec', spec='%s/A:%d@CA' % (SPEC, resnum)).get('centroid')
            if c:
                t = post('tug', spec='%s/A:%d' % (SPEC, resnum), target=[c[0] + 2.0, c[1], c[2]])
                show('tug A:%d' % resnum, 'n_tugged=%s' % t.get('n_tugged'))
                time.sleep(1)
                post('tug', spec='%s/A:%d' % (SPEC, resnum), release=True)
        # restraints
        rt = (inv('isolde restrain torsions', residues=SPEC).get('witness') or {})
        show('restrain torsions delta', rt.get('restraint_total_delta'))
        show('restraint_summary', post('restraint_summary', model=SPEC).get('restraints'))
        rel = (inv('isolde release torsions', residues=SPEC).get('witness') or {})
        show('release torsions enabled_delta', rel.get('restraint_enabled_delta'))
        # mdff toggle
        post('set_mdff', model=SPEC, enabled=False)
        show('set_mdff re-enable', [m['enabled'] for m in post('set_mdff', model=SPEC, enabled=True).get('mdff', [])])
        # checkpoint / revert
        inv('isolde sim', cmd='checkpoint'); show('checkpoint', 'ok')
        inv('isolde sim', cmd='revert'); show('revert', 'ok')
        inv('isolde sim', cmd='stop'); show('sim stop', 'requested')
        time.sleep(1)

    print('== a few more renders interleaved (no-op safety) ==')
    for _ in range(3):
        post('render', width=320, height=240)
    show('3x render', 'done (all no-op)')

    print('== final validation ==')
    r = inv('isolde validate ramachandran', model=SPEC).get('result', {})
    show('rama final', {k: r[k] for k in r if k.startswith('n_')})

    print('\nWORKOUT COMPLETE (%d failures). Now check the GUI: resize/drag/settle, dock/undock.' % nfail)


if __name__ == '__main__':
    main()
