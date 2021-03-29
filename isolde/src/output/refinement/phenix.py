


def write_phenix_restraints(file_name, 
        distance_restraints=None, 
        angle_restraints=None,
        dihedral_restraints=None):
    with open(file_name, 'wt') as outfile:
        outfile.write(
'''
geometry_restraints {
  edits {
'''
        )
        if distance_restraints is not None:
            for restraint in distance_restraints:
                write_distance_restraint(outfile, restraint)
        if angle_restraints is not None:
            for restraint in angle_restraints:
                write_angle_restraint(outfile, restraint)
        if dihedral_restraints is not None:
            for restraint in dihedral_restraints:
                write_dihedral_restraint(outfile, restraint)
        outfile.write(
'''
  }
}
'''            
        )
        

def phenix_selection_string(atom):
    cid = atom.residue.chain_id
    resnum = atom.residue.number
    name = atom.name
    altloc = atom.alt_loc
    out_str = 'chain {} and resid {} and name {}'.format(cid, resnum, name)
    if altloc != ' ':
        out_str += ' and altloc {}'.format(altloc)
    return out_str

def write_distance_restraint(outfile, restraint):
    atom1, atom2, distance, sigma, slack, top_out = restraint
    sp = ' '*4
    sp2 = ' '*6
    if slack is None:
        slack = 'None'
    else:
        slack = '{:.2f}'.format(slack)
    outfile.write(sp+'bond {\n')
    outfile.write(sp2+'action = *add\n')
    outfile.write(sp2+'atom_selection_1 = '+phenix_selection_string(atom1)+'\n')
    outfile.write(sp2+'atom_selection_2 = '+phenix_selection_string(atom2)+'\n')
    outfile.write(sp2+'symmetry_operation = None\n')
    outfile.write(sp2+'distance_ideal = {:.2f}\n'.format(distance))
    outfile.write(sp2+'sigma = {:.2f}\n'.format(sigma))
    outfile.write(sp2+'slack = {}\n'.format(slack))
    outfile.write(sp2+'limit = -1.0\n')
    outfile.write(sp2+'top_out = {}\n'.format(top_out))
    outfile.write(sp+'}\n')

def write_angle_restraint(outfile, restraint):
    #TODO
    pass

def write_dihedral_restraint(outfile, restraint):
    #TODO
    pass


_distances = {
    ('Mg','N'):     (2.2, 0.05),
    ('Fe','N'):     (2.05, 0.05),
    ('Fe','S'):     (2.30, 0.05)
}

_default_coordinating_atoms = {
    'HIS':  ('ND1','NE2'),
    'MET':  ('SD',)
}

def find_coordination_sites(model, residue_name, atom_name, coordinating=_default_coordinating_atoms, cutoff=3):
    from chimerax.geometry import find_closest_points
    import numpy
    residues = model.residues[model.residues.names==residue_name]
    atoms = residues.atoms[residues.atoms.names==atom_name]
    pairs = []
    for resname, atom_names in coordinating.items():
        crs = model.residues[model.residues.names==resname]
        catoms = crs.atoms[numpy.in1d(crs.atoms.names, atom_names)]
        i1, i2, near1 = find_closest_points(catoms.coords, atoms.coords, cutoff)
        for i, ii in zip(i1, near1):
            pairs.append([catoms[i], atoms[ii]])
    return pairs


def define_distance_restraints(pairs):
    distance_restraints = []
    for pair in pairs:
        names = tuple(sorted([a.element.name for a in pair]))
        dist, sigma = _distances[names]
        distance_restraints.append([*pair, dist, sigma, None, False])
    return distance_restraints
