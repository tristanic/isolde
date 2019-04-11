# @Author: Tristan Croll <tic20>
# @Date:   20-Dec-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 11-Apr-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



from math import radians
import numpy
_torsion_adjustments_chi1 = {
    'THR': radians(-120),
}
_torsion_adjustments_chi2 = {
    'TRP': radians(180),
}

def restrain_torsions_to_template(template_residues, restrained_residues,
    alignment_cutoff = 10, restrain_backbone=True, restrain_rotamers=True,
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
        * alignment_cutoff (default = 5):
            - distance cutoff (in Angstroms) for rigid-body alignment of model
              against  template. Residues with a CA RMSD greater than this
              value after alignment will not be restrained.
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
    # if len(template_residues) != len(restrained_residues):
    #     raise TypeError('Template and restrained residue arrays must be the same length!')
    template_us = template_residues.unique_structures
    if len(template_us) != 1:
        raise TypeError('Template residues must be from a single model!')
    template_model = template_us[0]
    restrained_us = restrained_residues.unique_structures
    if len(restrained_us) != 1:
        raise TypeError('Restrained residues must be from a single model!')
    restrained_model = restrained_us[0]

    if template_residues != restrained_residues:
        template_residues, restrained_residues = _get_template_alignment(
            template_residues, restrained_residues,
            cutoff_distance=alignment_cutoff,
            overlay_template=False
        )

    tdm = sx.get_proper_dihedral_mgr(template_model)
    rdrm = sx.get_proper_dihedral_restraint_mgr(restrained_model)
    names = ('phi','psi','omega','chi1','chi2','chi3','chi4')
    for name in names:
        for tr, rr in zip(template_residues, restrained_residues):
            td = tdm.get_dihedral(tr, name)
            rdr = rdrm.add_restraint_by_residue_and_name(rr, name)
            if td and rdr:
                # Due to naming conventions, some sidechain torsions need to be
                # rotated for best match to other residues
                if name == 'chi1':
                    target = (td.angle + _torsion_adjustments_chi1.get(tr.name, 0)
                        - _torsion_adjustments_chi1.get(rr.name, 0))
                elif name == 'chi2':
                    target = (td.angle + _torsion_adjustments_chi2.get(tr.name, 0)
                        - _torsion_adjustments_chi2.get(rr.name, 0))
                else:
                    target = td.angle
                rdr.target = target
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


def _get_template_alignment(template_residues, restrained_residues,
        cutoff_distance=5, overlay_template=False, always_raise_errors=False):
    ts = template_residues.unique_structures
    if len(ts) != 1:
            raise TypeError('Template residues must be from a single model! '
                'Residues are {} in {}'.format(template_residues.numbers, ','.join(s.id_string for s in template_residues.structures)))
    ts = ts[0]
    ms = restrained_residues.unique_structures
    if len(ms) != 1:
        raise TypeError('Residues to restrain must come from a single model!')
    ms = ms[0]
    same_model = False
    if ts == ms:
        # Matchmaker requires the two chains to be in separate models. Create a
        # temporary model here.
        same_model = True
        saved_ts = ts
        saved_trs = template_residues
        from chimerax.std_commands.split import molecule_from_atoms
        ts = molecule_from_atoms(ts, template_residues.atoms)
        template_residues = ts.residues




    from chimerax.match_maker.match import match, defaults
    result = match(ms.session, defaults['chain_pairing'], (ms, [ts]),
        defaults['matrix'], defaults['alignment_algorithm'],
        defaults['gap_open'], defaults['gap_extend'],
        cutoff_distance = cutoff_distance,
        domain_residues=(restrained_residues, template_residues),
        always_raise_errors=always_raise_errors)[0]
    if not overlay_template:
        # Return the template model to its original position
        ts.position = result[-1].inverse()*ts.position
    tr = result[0].residues
    rr = result[1].residues
    if same_model:
        # template_residues are in the temporary model. Need to find their
        # equivalents in the main model, and make sure they're in the right
        # order
        tr = saved_trs[template_residues.indices(tr)]
        ts.delete()
    return tr, rr

def restrain_atom_distances_to_template(template_residues, restrained_residues,
    protein=True, nucleic=True, custom_atom_names=[],
    distance_cutoff=8, alignment_cutoff=5, well_half_width = 0.05,
    kappa = 5, tolerance = 0.025, fall_off = 4):
    r'''
    Creates a "web" of adaptive distance restraints between nearby atoms,
    restraining one set of residues to the same spatial organisation as another.

    Args:
        * template_residues:
            - a list of :class:`chimerax.atomic.Residues` instances. If
              :attr:`restrained_residues` is not identical to
              :attr:`template_residues`, then each :class:`Residues` should be
              from a single chain. Residues need not be contiguous.
        * restrained_residues:
            - a list of :class:`chimerax.atomic.Residues` instances. Must be
              the same length as :attr:`template_residues`, with a 1:1
              correspondence between chains. Chains will be aligned individually
              to get the subset of matching residues, but original coordinates
              will be used for the purposes of assigning restraints.
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
        * alignment_cutoff (default = 5):
            - distance cutoff (in Angstroms) for rigid-body alignment of model
              against  template. Residues with a CA RMSD greater than this
              value after alignment will not be restrained.
        * well_half_width (default = 0.05):
            - distance range (as a fraction of the target distance) within which
              the restraint will behave like a normal harmonic restraint.
              The applied force will gradually taper off for any restraint
              deviating from (target + tolerance) by more than this amount.
        * kappa (default = 5):
            - defines the strength of each restraint when the current distance
              is within :attr:`well_half_width` of the target +/-
              :attr:`tolerance`. The effective spring constant is
              :math:`k=\frac{\kappa}{(\text{well\_half\_width}*\text{target distance})^2}`
              in :math:`kJ mol^{-1} nm^{-2}`.
        * tolerance (default = 0.025):
            - half-width (as a fraction of the target distance) of the "flat
              bottom" of the restraint profile. If
              :math:`abs(distance-target) < tolerance * target`,
              no restraining force will be applied.
        * fall_off (default = 6):
            - Sets the rate at which the energy function will fall off when the
              distance deviates strongly from the target, as a function of the
              target distance. The exponent on the energy term at large
              deviations from the target distance will be set as
              :math:`\alpha = -2 -\text{fall\_off} ln(\text{target})`. In other
              words, long-distance restraints are treated as less confident than
              short-distance ones.
    '''
    from chimerax.std_commands.align import IterationError
    from chimerax.isolde import session_extensions as sx
    if not protein and not nucleic and not len(custom_atom_names):
        raise TypeError('Nothing to restrain!')
    # if len(template_residues) != len(restrained_residues):
    #     raise TypeError('Template and restrained residue arrays must be the same length!')
    for trs in template_residues:
        template_us = trs.unique_structures
        if len(template_us) != 1:
            raise TypeError('Template residues must be from a single model! '
                'Residues are {} in {}'.format(trs.numbers, ','.join(s.id_string for s in trs.structures)))
    for rrs in restrained_residues:
        restrained_us = rrs.unique_structures
        if len(restrained_us) != 1:
            raise TypeError('Restrained residues must be from a single model!')
    restrained_model = restrained_us[0]
    log = restrained_model.session.logger

    adrm = sx.get_adaptive_distance_restraint_mgr(restrained_model)
    from chimerax.core.geometry import find_close_points, distance
    from chimerax.atomic import concatenate

    atom_names = []
    if protein:
        atom_names.extend(['CA', 'CB', 'CG', 'CG1', 'OG', 'OG1'])
    if nucleic:
        atom_names.extend(["OP1", "OP2", "C4'", "C2'", "O2", "O4", "N4", "N2", "O6", "N1", "N6"])
    atom_names.extend(custom_atom_names)

    import numpy
    atom_names = set(atom_names)

    def apply_restraints(trs, rrs):
        template_as = []
        restrained_as = []
        for tr, rr in zip(trs, rrs):
            ta_names = set(tr.atoms.names).intersection(atom_names)
            ra_names = set(rr.atoms.names).intersection(atom_names)
            common_names = list(ta_names.intersection(ra_names))
            template_as.extend([tr.find_atom(name) for name in common_names])
            restrained_as.extend([rr.find_atom(name) for name in common_names])
            # template_as.append(tr.atoms[numpy.in1d(tr.atoms.names, common_names)])
            # restrained_as.append(rr.atoms[numpy.in1d(rr.atoms.names, common_names)])
        from chimerax.atomic import Atoms
        template_as = Atoms(template_as)
        restrained_as = Atoms(restrained_as)

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
                dr.c = max(dist*well_half_width, 0.1)
                #dr.effective_spring_constant = spring_constant
                dr.kappa = kappa
                from math import log
                if dist < 1:
                    dr.alpha = -2
                else:
                    dr.alpha = -2 - fall_off * log(dist)
                dr.enabled = True


    if all(trs == rrs for trs, rrs in zip(template_residues, restrained_residues)):
        # If the template is identical to the model, we can just go ahead and restrain

        [apply_restraints(trs, rrs) for trs, rrs in zip(template_residues, restrained_residues)]
    else:
        # If we're not simply restraining the model to itself, we need to do an
        # alignment to get a matching pair of residue sequences
        while len(template_residues) != 0 and len(restrained_residues) != 0:
            found_trs = []
            found_rrs = []
            remain_trs = []
            remain_rrs = []
            for trs, rrs in zip(template_residues, restrained_residues):
                try:
                    ftrs, frrs = _get_template_alignment(
                        trs, rrs, cutoff_distance = alignment_cutoff,
                        overlay_template = False, always_raise_errors = True
                    )
                except IterationError:
                    log.info(("No sufficient alignment found for {} residues "
                        "in template chain {} and {} residues in restrained "
                        "model chain {}").format(len(trs), trs.chain_ids[0],
                        len(rrs), rrs.chain_ids[0]))
                    continue
                found_trs.append(ftrs)
                found_rrs.append(frrs)
                rtrs = trs.subtract(ftrs)
                rrrs = rrs.subtract(frrs)
                if len(rtrs) > 3 and len(rrrs) > 3:
                    remain_trs.append(rtrs)
                    remain_rrs.append(rrrs)
            template_residues = remain_trs
            restrained_residues = remain_rrs

            if len(found_trs):
                trs = concatenate(found_trs)
                rrs = concatenate(found_rrs)
                apply_restraints(trs, rrs)

def restrain_atom_pair_adaptive_distance(atom1, atom2, target, tolerance, kappa, c, alpha=-2):
    if not atom1.structure == atom2.structure:
        raise UserError('Both atoms must belong to the same model!')
    from chimerax.isolde import session_extensions as sx
    adrm = sx.get_adaptive_distance_restraint_mgr(atom1.structure)
    adr = adrm.add_restraint(atom1, atom2)
    adr.target = target
    adr.tolerance = tolerance
    adr.kappa = kappa
    adr.c = c
    adr.alpha = alpha
    adr.enabled = True
