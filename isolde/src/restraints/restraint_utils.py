from math import radians
import numpy
def restrain_torsions_to_template(template_residues, restrained_residues,
    restrain_backbone=True, restrain_rotamers=True,
    angle_cutoff=radians(10), spring_constant = 500):
    #from .. import session_extensions as sx
    from chimerax.isolde import session_extensions as sx
    if len(template_residues) != len(restrained_residues):
        raise TypeError('Template and restrained residue arrays must be the same length!')
    template_us = template_residues.unique_structures
    if len(template_us) != 1:
        raise TypeError('Template residues must be from a single model!')
    template_model = template_us[0]
    restrained_us = restrained_residues.unique_structures
    if len(restrained_us) != 1:
        raise TypeError('Restrained residues must be from a single model!')
    restrained_model = restrained_us[0]
    tdm = sx.get_proper_dihedral_mgr(template_model)
    rdrm = sx.get_proper_dihedral_restraint_mgr(restrained_model)
    names = ('phi','psi','omega','chi1','chi2','chi3','chi4')
    for name in names:
        for tr, rr in zip(template_residues, restrained_residues):
            td = tdm.get_dihedral(tr, name)
            rdr = rdrm.add_restraint_by_residue_and_name(rr, name)
            if td and rdr:
                rdr.target = td.angle
                rdr.spring_constant = spring_constant
                rdr.cutoff = angle_cutoff
                rdr.enabled = True

def restrain_ca_distances_to_template(template_residues, restrained_residues,
    distance_cutoff=8, spring_constant = 500):
    from chimerax.isolde import session_extensions as sx
    if len(template_residues) != len(restrained_residues):
        raise TypeError('Template and restrained residue arrays must be the same length!')
    template_us = template_residues.unique_structures
    if len(template_us) != 1:
        raise TypeError('Template residues must be from a single model!')
    restrained_us = restrained_residues.unique_structures
    if len(restrained_us) != 1:
        raise TypeError('Restrained residues must be from a single model!')
    restrained_model = restrained_us[0]
    template_cas = template_residues.atoms[template_residues.atoms.names=='CA']
    restrained_cas = restrained_residues.atoms[restrained_residues.atoms.names=='CA']
    template_coords = template_cas.coords
    drm = sx.get_distance_restraint_mgr(restrained_model)
    from chimerax.core.geometry import find_close_points, distance
    for i, rca1 in enumerate(restrained_cas):
        query_coord = numpy.array([template_coords[i]])
        indices = find_close_points(query_coord, template_coords, distance_cutoff)[1]
        indices = indices[indices !=i]
        for ind in indices:
            rca2 = restrained_cas[ind]
            dr = drm.add_restraint(rca1, rca2)
            dr.spring_constant = spring_constant
            dr.target = distance(query_coord[0], template_coords[ind])
            dr.enabled = True

def restrain_small_ligands(model, distance_cutoff=3.5, heavy_atom_limit=3, spring_constant=500,
    bond_to_carbon=False):
    '''
    Residues with a small number of heavy atoms can be problematic in MDFF if
    unrestrained, since if knocked out of density they tend to simply keep
    going. It is best to restrain them with distance restraints to suitable
    surrounding atoms or, failing that, to their starting positions.
    '''
    from chimerax.atomic import Residue, Residues
    residues = model.residues
    ligands = residues[residues.polymer_types==Residue.PT_NONE]
    small_ligands = Residues(
        [r for r in ligands if len(r.atoms[r.atoms.element_names!='H'])<heavy_atom_limit])
    from .. import session_extensions as sx
    drm = sx.get_distance_restraint_mgr(model)
    prm = sx.get_position_restraint_mgr(model)
    all_heavy_atoms = model.atoms[model.atoms.element_names!='H']
    if not bond_to_carbon:
        all_heavy_atoms = all_heavy_atoms[all_heavy_atoms.element_names != 'C']
    all_heavy_coords = all_heavy_atoms.coords

    from chimerax.core.geometry import find_close_points, distance
    for r in small_ligands:
        r_heavy_atoms = r.atoms[r.atoms.element_names != 'H']
        if not bond_to_carbon:
            r_non_carbon_atoms = r_heavy_atoms[r_heavy_atoms.element_names !='C']
            if not len(r_non_carbon_atoms):
                # No non-carbon heavy atoms. Apply position restraints
                prs = prm.get_restraints(r_heavy_atoms)
                prs.targets = prs.atoms.coords
                prs.spring_constants = spring_constant
                prs.enableds = True
                continue
            r_heavy_atoms = r_non_carbon_atoms
        r_indices = all_heavy_atoms.indices(r_heavy_atoms)
        r_coords = r_heavy_atoms.coords
        applied_drs = False
        for ra, ri, rc in zip(r_heavy_atoms, r_indices, r_coords):
            _, found_i = find_close_points([rc], all_heavy_coords, distance_cutoff)
            found_i = found_i[found_i!=ri]
            for fi in found_i:
                if fi in r_indices:
                    continue
                dr = drm.add_restraint(ra, all_heavy_atoms[fi])
                dr.spring_constant = spring_constant
                dr.target=distance(rc, all_heavy_coords[fi])
                dr.enabled = True
                applied_drs = True
        if not applied_drs:
            # Ligand is too far from the model for distance restraints. Use
            # position restraionts
            prs = prm.add_restraints(r_heavy_atoms)
            prs.targets = prs.atoms.coords
            prs.spring_constants = spring_constant
            prs.enableds = True
