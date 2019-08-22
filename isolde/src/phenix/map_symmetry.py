
def _parse_operator(lines):
    import numpy
    tf = numpy.zeros([3,4])
    done = False
    for i, line in enumerate(lines):
        if line.startswith('rota_matrix'):
            break
    else:
        raise RuntimeError('new_operator specified but no matrix found!')
    for j, line in enumerate(lines[i:i+3]):
        rot_operators = line.split()[1:4]
        tf[j,:3] = [float(o) for o in rot_operators]
    trn_operators = lines[i+3].split()[1:4]
    for j, op in enumerate(trn_operators):
        tf[j,3] = float(op)
    return (tf, i+4)



def parse_map_symmetry_file(filename):
    import numpy
    from collections import defaultdict
    ncs_map = defaultdict(list)
    with open(filename, 'rt') as infile:
        data = infile.read().split('\n')
    ncs_group_id = -1
    i = 0
    print(len(data))
    while i < len(data):
        line = data[i]
        if line.startswith('new_ncs_group'):
            ncs_group_id += 1
            i+=1
            continue
        if line.startswith('new_operator'):
            op, incr = _parse_operator(data[i+1:])
            ncs_map[ncs_group_id].append(op)
            i += (incr + 1)
            continue
        else:
            i+=1
            continue
    from chimerax.core.geometry import Places
    ret = {}
    for gid, ncslist in ncs_map.items():
        ret[gid] = Places(place_array=numpy.array(ncslist))
    return ret
