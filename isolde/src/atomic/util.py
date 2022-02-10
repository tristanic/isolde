rings = ('PHE','TYR') #,'TYS','PTR') <- non-standard amino acids currently throw a RuntimeError (25/10/2021)

def correct_pseudosymmetric_sidechain_atoms(session, residues):
    '''
    Protein sidechain atom names follow strict rules dictating the names of atoms 
    in symmetric sidechains (TYR, PHE, ASP, GLU) based on the preceding chi dihedral 
    angle. Additionally, MD forcefields make it very easy to "flip" an arginine sidechain
    by flipping just the NE atom - this renders the names of NH1 and NH2 reversed. Since we 
    don't want to limit mobility during model building, it's best to just automatically correct 
    these issues at suitable points. 
    '''
    session.logger.info('ISOLDE: Checking and correcting nomenclature for (pseudo)symmetric side chains...')
    import numpy
    from collections import defaultdict
    flipped = defaultdict(lambda: 0)
    for r in residues[numpy.in1d(residues.names, rings)]:
        if flip_if_necessary(r, ('CA','CB','CG',('CD1','CD2'))):
            flipped[r.structure] += 1
    for r in residues[residues.names=='ASP']:
        if not acid_is_substituted_or_protonated(r):
            if flip_if_necessary(r, ('CA','CB','CG',('OD1','OD2'))):
                flipped[r.structure] += 1
    for r in residues[residues.names=='GLU']:
        if not acid_is_substituted_or_protonated(r):
            if flip_if_necessary(r, ('CB','CG','CD',('OE1','OE2'))):
                flipped[r.structure] += 1
    for r in residues[residues.names=='ARG']:
        if flip_if_necessary(r, ('CD','NE','CZ',('NH2','NH1'))):
            flipped[r.structure] += 1
    for m, count in flipped.items():
        session.logger.info(f'ISOLDE: Corrected atom nomenclature of {count} residues in model #{m.id_string} to IUPAC-IUB standards.')
        

def any_atom_restrained(model, atoms):
    import numpy
    from chimerax.isolde import session_extensions as sx
    prm = sx.get_position_restraint_mgr(model, create=False)
    drm = sx.get_distance_restraint_mgr(model, create=False)
    adrm = sx.get_adaptive_distance_restraint_mgr(model, create=False)
    if prm is not None:
        prs = prm.get_restraints(atoms)
        if numpy.any(prs.enableds):
            return True
    if drm is not None:
        drs = drm.atoms_restraints(atoms)
        if numpy.any(drs.enableds):
            return True
    if adrm is not None:
        adrs = adrm.atoms_restraints(atoms)
        if numpy.any(adrs.enableds):
            return True
    return False


    
def flip_if_necessary(residue, chi_atom_names):
    if residue.is_missing_heavy_template_atoms():
        return
    chis = [[*chi_atom_names[:3], chi_atom_names[3][i]] for i in range(2)]
    dihedrals = []
    from chimerax.geometry import dihedral
    for chi_atom_names in chis:
        atoms = [residue.find_atom(name) for name in chi_atom_names]
        dihedrals.append(dihedral(*[a.coord for a in atoms]))
    if abs(dihedrals[1]) < abs(dihedrals[0]):
        if swap_equivalent_atoms(residue):
            correct_dihedral_restraint(residue)
            correct_position_and_distance_restraints(residue)
            return True
    return False

def acid_is_substituted_or_protonated(residue):
    if residue.name == 'ASP':
        pos = 'D'
    elif residue.name == 'GLU':
        pos = 'E'
    else:
        raise TypeError('This method is only applicable to ASP or GLU residues!')
    for i in (1,2):
        o = residue.find_atom(f'O{pos}{i}')
        if o is None:
            return True
        if len(o.neighbors) > 1:
            return True
    return False

equivalent_heavy_atoms = {
    'aromatic': {'CD1':'CD2', 'CE1':'CE2'},
    'ASP':  {'OD1': 'OD2'},
    'GLU':  {'OE1': 'OE2'},
    'ARG':  {'NH1': 'NH2'}
}

def swap_equivalent_atoms(residue):
    rname = residue.name
    if rname in rings:
        rname = 'aromatic'
    equivalent_pairs = equivalent_heavy_atoms[rname]
    from chimerax.atomic import Atoms
    paired_atoms = []
    for a1name, a2name in equivalent_pairs.items():
        a1, a2 = [residue.find_atom(name) for name in (a1name, a2name)]
        if a1 is None or a2 is None:
            # Bail out if any equivalent atom is missing
            return False
        pair = []
        for a in (a1, a2):
            group = [a]
            for n in a.neighbors:
                if n.residue != a.residue:
                    # Don't flip if any equivalent atom is bonded to something else
                    return False
                if n.element.name == 'H':
                    group.append(n)
            pair.append(Atoms(group))
        paired_atoms.append(pair)
    for pair in paired_atoms:
        pair[0].coords, pair[1].coords = pair[1].coords, pair[0].coords
    return True


chi_names = {
    'aromatic': 'chi2',
    'ASP': 'chi2',
    'GLU': 'chi3',
    'ARG': None
}

def correct_dihedral_restraint(residue):
    from chimerax.isolde import session_extensions as sx
    from math import radians
    cutoff = radians(90)
    model = residue.structure
    from math import pi
    rname = residue.name
    if rname in rings:
        rname='aromatic'
    chi_name = chi_names.get(rname, None)
    if chi_name is None:
        return
    pdrm = sx.get_proper_dihedral_restraint_mgr(model, create=False)
    apdrm = sx.get_adaptive_dihedral_restraint_mgr(model, create=False)
    if pdrm is not None:
        pdr = pdrm.get_restraint_by_residue_and_name(residue, chi_name)
        if pdr is not None:
            # Only correct if it puts the current angle closer to the target
            if abs(pdr.offset) > cutoff:
                pdr.target += pi
    if apdrm is not None:
        apdr = apdrm.get_restraint_by_residue_and_name(residue, chi_name)
        if apdr is not None:
            if abs(apdr.offset) > cutoff:
                apdr.target += pi

    

def correct_position_and_distance_restraints(residue):
    rname = residue.name
    model = residue.structure
    from chimerax.isolde import session_extensions as sx
    if rname in rings:
        rname = 'aromatic'
    pairs = equivalent_heavy_atoms[rname]
    apairs = {residue.find_atom(a1name): residue.find_atom(a2name) for a1name, a2name in pairs.items()}
    prm = sx.get_position_restraint_mgr(model, create=False)
    if prm is not None:
        for a1, a2 in apairs.items():
            p1, p2 = [prm.get_restraint(a) for a in [a1, a2]]
            params = {}
            for p in [p1, p2]:
                if p is not None:
                    params[p] = [p.enabled, p.spring_constant, p.target]
            if not len(params):
                continue
            if p1 is None and p2 is not None:
                p1 = prm.add_restraint(a1)
            elif p2 is None and p1 is not None:
                p2 = prm.add_restraint(a2)
            p2params = params.get(p2, None)
            if p2params:
                p1.enabled, p1.spring_constant, p1.target = p2params
            p1params = params.get(p1, None)
            if p1params:
                p2.enabled, p2.spring_constant, p2.target = p1params

    drm = sx.get_distance_restraint_mgr(model, create=False)
    adrm = sx.get_adaptive_distance_restraint_mgr(model, create=False)
    if drm is not None:
        restraint_params = {}
        for a1, a2 in apairs.items():
            for a in (a1, a2):
                restraints = drm.atom_restraints(a)
                restraint_params[a] = [[r.atoms, r.enabled, r.target, r.spring_constant] for r in restraints]
                restraints.enableds=False
            for rparams in restraint_params[a1]:
                atoms = rparams[0]
                other = [a for a in atoms if a!= a1][0]
                if other != a2:
                    atoms = [a2, other]
                restraint = drm.add_restraint(*atoms)
                restraint.enabled, restraint.target, restraint.spring_constant = rparams[1:]
            for rparams in restraint_params[a2]:
                atoms = rparams[0]
                other = [a for a in atoms if a!= a2][0]
                if other != a1:
                    atoms = [a1, other]
                restraint = drm.add_restraint(*atoms)
                restraint.enabled, restraint.target, restraint.spring_constant = rparams[1:]
    
    if adrm is not None:
        restraint_params = {}
        for a1, a2 in apairs.items():
            for a in (a1, a2):
                restraints = adrm.atom_restraints(a)
                restraint_params[a] = [[r.atoms, r.enabled, r.target, r.kappa, r.alpha, r.tolerance, r.c] for r in restraints]
                restraints.enableds=False
            for rparams in restraint_params[a1]:
                atoms = rparams[0]
                other = [a for a in atoms if a!= a1][0]
                if other != a2:
                    atoms = [a2, other]
                r = adrm.add_restraint(*atoms)
                r.enabled, r.target, r.kappa, r.alpha, r.tolerance, r.c = rparams[1:]
            for rparams in restraint_params[a2]:
                atoms = rparams[0]
                other = [a for a in atoms if a!= a2][0]
                if other != a1:
                    atoms = [a1, other]
                r = adrm.add_restraint(*atoms)
                r.enabled, r.target, r.kappa, r.alpha, r.tolerance, r.c = rparams[1:]
        





