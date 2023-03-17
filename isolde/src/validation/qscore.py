def unit_sphere_vertices(num_vertices):
    import numpy
    from chimerax.surface.shapes import sphere_geometry2
    sphere_vertices = sphere_geometry2(2*num_vertices-4)[0] # 128 points evenly distributed around a unit sphere centred on (0,0,0)
    return sphere_vertices

def min_max_d(v):
    m = v.data.full_matrix()
    mean, sd, _ = v.mean_sd_rms()
    max_d = min (mean+sd*10, m.max())
    min_d = max (mean-sd, m.min())
    return min_d, max_d

_numpy_print_options = {
    'precision':3, 
    'suppress':True
    }


def q_score(residues, volume, ref_sigma=0.6, points_per_shell = 8, max_rad = 2.0, step=0.1, num_test_points=128, clustering_iterations=5, include_h = False, debug = False, draw_points=False):
    from chimerax.geometry import find_close_points, find_closest_points, Places
    import numpy
    from math import floor
    if len(residues.unique_structures) != 1:
        from chimerax.core.errors import UserError
        raise UserError('All residues must be from the same model!')
    session = residues.unique_structures[0].session
    from .._kmeans import spherical_k_means
    matrix, xform = volume.matrix_and_transform(None, 'all', (1,1,1))
    from chimerax.map_data import interpolate_volume_data
    min_d, max_d = min_max_d(volume)
    sphere8_vertices = unit_sphere_vertices(8)
    ref_sphere_vertices = unit_sphere_vertices(num_test_points)

    a = max_d - min_d
    b = min_d

    radii = numpy.arange(0, max_rad+step/2, step)
    reference_gaussian = a * numpy.exp( -0.5 * (radii/ref_sigma)**2) + b

    query_atoms = residues.atoms
    all_atoms = residues.unique_structures[0].atoms
    if not include_h:
        query_atoms = query_atoms[query_atoms.element_names != 'H']
        all_atoms = all_atoms[all_atoms.element_names != 'H']



    if draw_points:
        from chimerax.core.models import Model, Drawing
        from chimerax.surface.shapes import sphere_geometry2
        from chimerax.core.colors import random_colors
        d = Drawing('shell points')
        v, n, t = sphere_geometry2(24)
        v *= (step*0.25)
        d.set_geometry(v, n, t)
        dm = Model('shell points', session)
        dm.add_drawing(d)
        positions = []
        base_colors = random_colors(len(query_atoms))
        colors = []

    all_coords = all_atoms.scene_coords
    query_coords = query_atoms.scene_coords

    r0_vals = volume.interpolated_values(query_coords)
    
    num_shells = int(floor(max_rad/step))


    q_scores = []
    for i,a in enumerate(query_atoms):
        if draw_points:
            color = base_colors[i]
        a_coord = a.scene_coord
        _, nearby_i = find_close_points([a_coord], all_coords, max_rad*3)
        nearby_a = all_atoms[nearby_i]
        ai = nearby_a.index(a)
        nearby_coords = nearby_a.scene_coords
        shell_rad = step
        local_d_vals = {}
        j = 1
        while shell_rad < max_rad+step/2:
            if shell_rad < 0.7: # about half a C-C bond length
                # Try the quick way first (should succeed for almost all cases unless geometry is seriously wonky)
                local_8 = (sphere8_vertices*shell_rad) + a_coord
                i1, i2, near1 = find_closest_points(local_8, nearby_coords, shell_rad*1.5)
                closest = near1
                candidates = i1[closest==ai]
                if len(candidates)==8:
                    d_vals = interpolate_volume_data(local_8, xform, matrix, 'linear')[0]
                    if draw_points:
                        positions.append(local_8)
                        colors.append(numpy.array([color]*8))
                    local_d_vals[j] = d_vals
                    shell_rad += step
                    j+=1
                    continue
            
            local_sphere = (ref_sphere_vertices*shell_rad) + a_coord
            i1, i2, near1 = find_closest_points(local_sphere, nearby_coords, shell_rad*1.5)
            closest = near1
            candidates = i1[closest==ai]
            if len(candidates) == 0:
                if debug:
                    id_string = f'#{a.residue.structure.id_string}/{a.residue.chain_id}:{a.residue.number}@{a.name}'
                    print(f'WARNING: no shell points found for atom {id_string} at radius {shell_rad:.1f}!')
            elif len(candidates) < points_per_shell:
                points = local_sphere[candidates]
                d_vals = interpolate_volume_data(points, xform, matrix, 'linear')[0]
                local_d_vals[j] = d_vals
            else:
                points = local_sphere[candidates]
                labels, closest = spherical_k_means(points, a_coord, points_per_shell, clustering_iterations)
                if debug:
                    with numpy.printoptions(**_numpy_print_options):
                        print(f'Points: {points}')
                        print(f'Labels: {labels}')
                        print(f'Closest indices: {closest}')
                points = points[closest]
                d_vals = interpolate_volume_data(points, xform, matrix, 'linear')[0]
                local_d_vals[j] = d_vals
            if draw_points:
                if len(candidates) != 0:
                    positions.append(points)
                    colors.append(numpy.array([color]*len(points)))
            shell_rad += step
            j+=1
        from chimerax.map_fit import overlap_and_correlation
        measured = numpy.concatenate([[r0_vals[i]]*8, *local_d_vals.values()])
        ref = numpy.concatenate([[reference_gaussian[0]]*8, *[[reference_gaussian[j]] * len(vals) for j,vals in local_d_vals.items()]])
        q = overlap_and_correlation(ref, measured)[2]
        if debug:
            with numpy.printoptions(**_numpy_print_options):
                print(f'Measured: {measured}')
                print(f'Ref: {ref}')
                print(f'q: {q:.3f}')
        q_scores.append(q)      

    if draw_points:
        positions = numpy.concatenate(positions)
        colors = numpy.concatenate(colors)
        shift_and_scale = numpy.ones((len(positions),4), dtype=numpy.float32)
        shift_and_scale[:,:3] = positions
        positions = Places(shift_and_scale=shift_and_scale)
        d.positions = positions
        d.colors = colors
        session.models.add([dm])
        
    q_scores = numpy.array(q_scores)
    residue_scores = {}
    for r in residues:
        indices = query_atoms.indices(r.atoms)
        indices = indices[indices!=-1]
        if len(indices) > 0:
            rscores = q_scores[indices]
            residue_scores[r] = (rscores.mean(), rscores.min())


    return residue_scores, q_scores


def test_q_score(session):
    from chimerax.core.commands import run
    m = run(session, 'open 6out')[0]
    v = run(session, 'open 20205 from emdb')[0]
    residue_map, atom_scores = q_score(m.residues, v)
    from matplotlib import pyplot as plt
    fig = plt.figure()
    plots = fig.subplots(len(m.chains))
    for i,c in enumerate(m.chains):
        plot = plots[i]
        residues = c.residues
        x = []
        y = []
        x = [r.number for r in residues if r is not None]
        y = [residue_map[r][0] for r in residues if r is not None]
        plot.plot(x, y)
        plot.set_xlabel('Residue number')
        plot.set_ylabel('Q Score')
        plot.set_title(f'Chain {c.chain_id}')
    fig.tight_layout()
    fig.show()












    





