# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



def pick_closest_to_line(session, mx, my, atoms, cutoff, displayed_only = True, hydrogens = False):
    '''
    Pick the atom coming closest to the ray projected from the mouse pointer
    away from the camera. Only atoms found between the near and far clipping
    planes and within cutoff of the line will be considered. Optionally the
    selection can be further limited to include only displayed atoms and/or
    exclude hydrogens.
    '''
    closest = None
    if atoms is None:
        return None
    xyz1, xyz2 = session.main_view.clip_plane_points(mx, my)
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
    if displayed_only:
        atoms = atoms.filter(atoms.visibles)
    if not hydrogens:
        atoms = atoms.filter(atoms.element_names != 'H')
    atomic_coords = atoms.scene_coords
    from chimerax.geometry import find_close_points
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

def pick_closest_to_point(session, xyz, atoms, cutoff, displayed_only = True, hydrogens = False):
    '''
    Pick the atom closest to an (x,y,z) point in scene coordinates. Optionally the
    selection can be limited to include only displayed atoms and/or
    exclude hydrogens.
    '''
    closest = None
    import numpy
    if displayed_only:
        atoms = atoms.filter(atoms.displays)
    if not hydrogens:
        atoms = atoms.filter(atoms.element_names != 'H')
    atomic_coords = atoms.scene_coords

    from chimerax.geometry import find_closest_points
    ignore1, ignore2, nearest = find_closest_points([xyz], atomic_coords, cutoff)
    if nearest is not None:
        if len(nearest):
            return atoms[nearest[0]]
    return
