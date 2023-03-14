def unit_sphere_vertices(num_vertices):
    import numpy
    from chimerax.surface.shapes import sphere_geometry2
    sphere_vertices = sphere_geometry2(2*num_vertices-4)[0] # 128 points evenly distributed around a unit sphere centred on (0,0,0)
    #sphere8 =  1.1547 * (numpy.array([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],dtype=numpy.float32)-0.5) # 8 points evenly distributed around a unit sphere centred on (0,0,0)
    return sphere_vertices

def min_max_d(v):
    m = v.data.full_matrix()
    mean, sd, _ = v.mean_sd_rms()
    max_d = min (mean+sd*10, m.max())
    min_d = max (mean-sd, m.min())
    return min_d, max_d

def q_score(residues, volume, ref_sigma=0.6, points_per_shell = 8, max_rad = 2.0, step=0.1, num_test_points=128, include_h = False):
    from chimerax.geometry import find_close_points, find_closest_points
    from chimerax.atomic import concatenate
    import numpy
    from math import floor
    from .._kmeans import spherical_k_means
    matrix, xform = volume.matrix_and_transform(None, 'all', (1,1,1))
    from chimerax.map_data import interpolate_volume_data
    min_d, max_d = min_max_d(volume)
    ref_sphere_vertices = unit_sphere_vertices(num_test_points)

    a = max_d - min_d
    b = min_d

    radii = numpy.arange(0, max_rad+step/2, step)
    reference_gaussian = a * numpy.exp( -0.5 * (radii/ref_sigma)**2) + b

    all_atoms = residues.atoms
    if not include_h:
        all_atoms = all_atoms[all_atoms.element_names != 'H']
    
    atom_q_scores = numpy.zeros(len(all_atoms), dtype=numpy.float32)
    residue_q_scores = numpy.zeros(len(residues), dtype=numpy.float32)

    all_coords = all_atoms.scene_coords

    r0_vals = volume.interpolated_values(all_coords)
    
    num_shells = int(floor(max_rad/step))
    r_vals = numpy.empty((len(all_atoms), num_shells+1))
    r_vals[:,0] = r0_vals

    for i,a in enumerate(all_atoms):
        a_coord = a.scene_coord
        _, nearby_i = find_close_points([a_coord], all_coords, max_rad*2+.1)
        nearby_a = all_atoms[nearby_i]
        ai = nearby_a.index(a)
        # nearby_a = concatenate([nearby_a[:ai],nearby_a[ai+1:]])
        nearby_coords = nearby_a.scene_coords
        shell_rad = step
        j = 1
        while shell_rad < max_rad+step/2:
            local_sphere = (ref_sphere_vertices*shell_rad) + a_coord
            i1, i2, near1 = find_closest_points(local_sphere, nearby_coords, shell_rad)
            closest = near1
            candidates = closest[closest==ai]
            if len(candidates) == 0:
                # id_string = f'#{a.residue.structure.id_string}/{a.residue.chain_id}:{a.residue.number}@{a.name}'
                # print(f'WARNING: no shell points found for atom {id_string} at radius {shell_rad:.1f}!')
                r_vals[i,j] = 0
            elif len(candidates) < points_per_shell:
                points = local_sphere[candidates]
                r_vals[i,j] = interpolate_volume_data(points, xform, matrix, 'linear')[0].mean()
            else:
                points = local_sphere[candidates]
                labels, closest = spherical_k_means(points, points_per_shell, 5, True)
                points = points[closest]
                r_vals[i,j] = interpolate_volume_data(points, xform, matrix, 'linear')[0].mean()
            shell_rad += step
            j+=1
    

    from chimerax.map_fit import overlap_and_correlation
    q_scores = numpy.array([overlap_and_correlation(reference_gaussian, r_vals[i,:])[2] for i in range(r_vals.shape[0])])
    residue_scores = {}
    for r in residues:
        indices = all_atoms.indices(r.atoms)
        indices = indices[indices!=-1]
        if len(indices) > 0:
            rscores = q_scores[indices]
            residue_scores[r] = (rscores.mean(), rscores.min())

    #q_scores = (r_vars*rg_var)/(r_norms*rg_norm)


    return residue_scores


def find_best_eight_points(points, sphere8, shell_radius):
    '''
    Given a set of points distributed on a sphere, try to find the eight most widely 
    spaced.
    '''
    from chimerax.geometry import find_closest_points
    import numpy

    i1, i2, near1 = find_closest_points(sphere8, points, shell_radius)
    if len(numpy.unique(near1)):
        return points[i2][near1]











    





