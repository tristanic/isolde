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
