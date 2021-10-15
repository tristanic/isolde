def correct_pseudosymmetric_sidechain_atoms(session, residues):
    '''
    Protein sidechain atom names follow strict rules dictating the names of atoms 
    in symmetric sidechains (TYR, PHE, ASP, GLU) based on the preceding chi dihedral 
    angle. Additionally, MD forcefields make it very easy to "flip" an arginine sidechain
    by flipping just the NE atom - this renders the names of NH1 and NH2 reversed. Since we 
    don't want to limit mobility during model building, it's best to just automatically correct 
    these issues at suitable points. 
    '''
    session.logger.info('Checking and correcting nomenclature for (pseudo)symmetric side chains...')
    rings = ('PHE','TYR','TYS','PTR')
    import numpy
    for r in residues[numpy.in1d(residues.names, rings)]:
        flip_if_necessary(session, r, ('CA','CB','CG','CD1'))
    for r in residues[residues.names=='ASP']:
        flip_if_necessary(session, r, ('CA','CB','CG','OD1'))
    for r in residues[residues.names=='GLU']:
        flip_if_necessary(session, r, ('CB','CG','CD1','OE1'))
    for r in residues[residues.names=='ARG']:
        flip_if_necessary(session, r, ('CD','NE','CZ','NH2'))
        

    
def flip_if_necessary(session, residue, chi_atom_names):
    atoms = [residue.find_atom(name) for name in chi_atom_names]
    if not all(atoms):
        # residue is incomplete
        return
    from chimerax.geometry import dihedral
    angle = dihedral(*[a.coord for a in atoms])
    if abs(angle) > 90:
        flip_sidechain_on_bond(session, atoms[1:3])



def flip_sidechain_on_bond(session, atoms):
    a, b = atoms
    for bond in a.bonds:
        if b in bond.atoms:
            break
    else:
        raise TypeError('Atoms must be bonded!')
    moving_atom = bond.smaller_side
    fixed_atom = bond.other_atom(moving_atom)
    m, f = moving_atom.coord, fixed_atom.coord
    moving_atoms = bond.side_atoms(moving_atom)
    if len(moving_atoms.unique_residues) != 1:
        # sidechain is bonded to another residue. Bail out.
        return
    residue = moving_atom.residue
    session.logger.info(f'Correcting sidechain atoms of {residue.name} #{residue.structure.id_string}/{residue.chain_id}:{residue.number}{residue.insertion_code} to IUPAC-IUB standards.')
    from chimerax.geometry import z_align, rotation
    za = z_align(m, f)
    tf = za.inverse() * rotation((0,0,-1), 180) * za
    moving_atoms.transform(tf)







