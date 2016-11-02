def pick_closest_to_line(session, mx, my, atoms, cutoff):
    closest = None
    xyz1, xyz2 = session.view.clip_plane_points(mx, my)
    import numpy
    # Create an array of coordinates with spacing cutoff/2
    length = numpy.linalg.norm(xyz2-xyz1)
    numpoints = numpy.ceil(length/cutoff*2).astype(int)
    xvals = numpy.linspace(xyz1[0],xyz2[0],num=numpoints)
    yvals = numpy.linspace(xyz1[1],xyz2[1],num=numpoints)
    zvals = numpy.linspace(xyz1[2],xyz2[2],num=numpoints)
    xyzlist = []
    for xyz in zip(xvals, yvals, zvals):
        xyzlist.append(xyz)
    xyzlist = numpy.array(xyzlist)
    atomic_coords = atoms.coords
    from chimerax.core.geometry import find_close_points
    line_indices, atom_indices = find_close_points(xyzlist, atomic_coords, cutoff)
    line_shortlist = xyzlist[line_indices]
    ac_shortlist = atomic_coords[atom_indices]
    atom_shortlist = atoms[atom_indices]
    min_dist = cutoff
    for lxyz in line_shortlist:
        for axyz, atom in zip(ac_shortlist, atom_shortlist):
            d = numpy.linalg.norm(axyz-lxyz)
            if d < min_dist:
                closest = atom
                min_dist = d
    return closest
            
    
