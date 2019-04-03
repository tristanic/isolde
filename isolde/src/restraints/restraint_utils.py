# @Author: Tristan Croll <tic20>
# @Date:   20-Dec-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 03-Apr-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



from math import radians
import numpy
def restrain_torsions_to_template(template_residues, restrained_residues,
    restrain_backbone=True, restrain_rotamers=True,
    angle_cutoff=radians(10), spring_constant = 500):
    '''
    (EXPERIMENTAL)

    Restrain all phi, psi, omega and/or chi dihedrals in `restrained_residues`
    to  match their counterparts (if present) in template_residues. There must
    be a 1:1 correspondence between the residues in the two arrays. Note that
    this algorithm is still quite simplistic and doesn't yet attempt to handle
    "gotchas" when the restrained and template residues are different (e.g.
    restraining a VAL to be isosteric with a THR requires rotating the chi1
    torsion by 120 degrees).

    Args:
        * template_residues:
            - a :class:`chimerax.atomic.Residues` instance. All residues must be
              from a single model, but need no be contiguous
        * restrained_residues:
            - a :class:`chimerax.atomic.Residues` instance. All residues must be
              from a single model (which may or may not be the same model as for
              `template_residues`). May be the same array as `template_residues`
              (which will just restrain all torsions to their current angles).
        * restrain_backbone (default = `True`):
            - if `True`, all phi, psi and omega dihedrals in
              `restrained_residues` that  exist in both `restrained_residues`
              and `template_residues` will be restrained to the angles in
              `template_residues`
        * restrain_rotamers (default = `True`):
            - if `True`, all chi dihedrals in `restrained_residues` that  exist
              in both `restrained_residues` and `template_residues` will be
              restrained to the angles in `template_residues`
        * angle_cutoff (default = pi/18 (10 degrees)):
            - the deviation from the target angle in radians below which no
              restraining force will be applied to a given torsion.
        * spring_constant (default = 500):
            - strength of each restraint, in :math:`kJ mol^{-1} rad^{-2}`
    '''
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
    '''
    Creates a "web" of distance restraints between nearby CA atoms, restraining
    one set of residues to the same spatial organisation as another.

    Args:
        * template_residues:
            - a :class:`chimerax.atomic.Residues` instance. All residues must be
              from a single model, but need no be contiguous
        * restrained_residues:
            - a :class:`chimerax.atomic.Residues` instance. All residues must be
              from a single model (which may or may not be the same model as for
              `template_residues`). May be the same array as `template_residues`
              (which will just restrain all distances to their current values).
        * distance_cutoff (default = 8):
            - for each CA atom in `restrained_residues`, a distance restraint
              will be created between it and every other CA atom where the
              equivalent atom in `template_residues` is within `distance_cutoff`
              of its template equivalent.
        * spring_constant (default = 500):
            - the strength of each restraint, in :math:`kJ mol^{-1} nm^{-2}`
    '''
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

    Args:
        * model:
            - a :class:`chimerax.atomic.AtomicStructure` instance
        * distance_cutoff (default = 3.5):
            - radius in Angstroms to look for candidate heavy atoms for distance
              restraints. If no candidates are found, a position restraint will
              be applied instead.
        * heavy_atom_limit (default = 3):
            - Only residues with a number of heavy atoms less than or equal to
              `heavy_atom_limit` will be restrained
        * spring_constant (default = 500):
            - strength of each restraint, in :math:`kJ mol^{-1} nm^{-2}`
        * bond_to_carbon (default = `False`):
            - if `True`, only non-carbon heavy atoms will be restrained using
              distance restraints.
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


def restrain_atom_distances_to_template(template_residues, restrained_residues,
    protein=True, nucleic=True, custom_atom_names=[],
    distance_cutoff=8, restrained_range = 0.2,
    spring_constant = 50, tolerance = 0.05, fall_off = 2):
    r'''
    Creates a "web" of adaptive distance restraints between nearby atoms,
    restraining one set of residues to the same spatial organisation as another.

    Args:
        * template_residues:
            - a :class:`chimerax.atomic.Residues` instance. All residues must be
              from a single model, but need no be contiguous
        * restrained_residues:
            - a :class:`chimerax.atomic.Residues` instance. All residues must be
              from a single model (which may or may not be the same model as for
              `template_residues`). May be the same array as `template_residues`
              (which will just restrain all distances to their current values).
        * protein (default = True):
            - Restrain protein conformation? If True, a pre-defined set of
              useful "control" atoms (CA plus the first two heavy atoms along
              each sidechain) will be added to the restraint network.
        * nucleic (default = True):
            - Restrain nucleic acid conformation? If True, key atoms defining
              the nucleic acid backbone and base pairing will be added to the
              restraint network.
        * custom_atom_names(default = empty list):
            - Provide the names of any other atoms you wish to restrain (e.g.
              ligand atoms) here.
        * distance_cutoff (default = 8):
            - for each CA atom in `restrained_residues`, a distance restraint
              will be created between it and every other CA atom where the
              equivalent atom in `template_residues` is within `distance_cutoff`
              of its template equivalent.
        * restrained_range (default = 0.2):
            - distance range (as a fraction of the target distance) within which
              the restraint will behave like a normal harmonic restraint.
              The applied force will gradually taper off for any restraint
              deviating from (target + tolerance) by more than this amount.
        * spring_constant (default = 50):
            - the strength of each restraint when the current distance is
              within :attr:`restrained_range` of the target +/-
              :attr:`tolerance`, in :math:`kJ mol^{-1} nm^{-2}`.
        * tolerance (default = 0.05):
            - half-width (as a fraction of the target distance) of the "flat
              bottom" of the restraint profile. If
              :math:`abs(distance-target) < tolerance * target`,
              no restraining force will be applied.
        * fall_off (default = 2):
            - Sets the rate at which the energy function will fall off when the
              distance deviates strongly from the target, as a function of the
              target distance. The exponent on the energy term at large
              deviations from the target distance will be set as
              :math:`\alpha = -2 -\text{fall\_off} ln(\text{target})`. In other
              words, long-distance restraints are treated as less confident than
              short-distance ones.
    '''
    from chimerax.isolde import session_extensions as sx
    if not protein and not nucleic and not len(custom_atom_names):
        raise TypeError('Nothing to restrain!')
    if len(template_residues) != len(restrained_residues):
        raise TypeError('Template and restrained residue arrays must be the same length!')
    template_us = template_residues.unique_structures
    if len(template_us) != 1:
        raise TypeError('Template residues must be from a single model!')
    restrained_us = restrained_residues.unique_structures
    if len(restrained_us) != 1:
        raise TypeError('Restrained residues must be from a single model!')
    restrained_model = restrained_us[0]
    adrm = sx.get_adaptive_distance_restraint_mgr(restrained_model)
    from chimerax.core.geometry import find_close_points, distance

    atom_names = []
    if protein:
        atom_names.extend(['CA', 'CB', 'CD', 'CD1', 'CG1', 'OG'])
    if nucleic:
        atom_names.extend(["OP1", "OP2", "C4'", "C2'", "O2", "O4", "N4", "N2", "O6", "N1", "N6"])
    atom_names.extend(custom_atom_names)

    import numpy
    atom_names = numpy.array(atom_names)

    template_as = template_residues.atoms[numpy.in1d(template_residues.atoms.names, atom_names)]
    restrained_as = restrained_residues.atoms[numpy.in1d(restrained_residues.atoms.names, atom_names)]
    template_coords = template_as.coords
    for i, ra1 in enumerate(restrained_as):
        query_coord = numpy.array([template_coords[i]])
        indices = find_close_points(query_coord, template_coords, distance_cutoff)[1]
        indices = indices[indices !=i]
        for ind in indices:
            ra2 = restrained_as[ind]
            if ra1.residue == ra2.residue:
                continue
            try:
                dr = adrm.add_restraint(ra1, ra2)
            except ValueError:
                continue
            dist = distance(query_coord[0], template_coords[ind])
            dr.tolerance = tolerance * dist
            dr.target = dist
            dr.c = max(dist*restrained_range, 0.25)
            dr.effective_spring_constant = spring_constant
            from math import log
            dr.alpha = -2 - fall_off * log(dist)
            dr.enabled = True
