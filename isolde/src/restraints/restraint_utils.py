# @Author: Tristan Croll <tic20>
# @Date:   20-Dec-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 10-Dec-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



from math import radians
import numpy
_torsion_adjustments_chi1 = {
    'THR': radians(-120),
}
_torsion_adjustments_chi2 = {
    'TRP': radians(180),
}

KAPPA_MAX=30

from chimerax.core.errors import UserError

def restrain_torsions_to_template(session, template_residues, restrained_residues,
    restrain_backbone=True, restrain_sidechains=True,
    kappa=10, alpha=0.2, spring_constant=100, identical_sidechains_only=True,
    adjust_for_confidence=False, confidence_type='plddt'):
    r'''
    Restrain all phi, psi, omega and/or chi dihedrals in `restrained_residues`
    to  match their counterparts (if present) in template_residues.

    Args:
        * template_residues:
            - a :class:`chimerax.atomic.Residues` instance. All residues must be
              from a single model, but need no be contiguous
        * restrained_residues:
            - a :class:`chimerax.atomic.Residues` instance. All residues must be
              from a single model (which may or may not be the same model as for
              `template_residues`). May be the same array as `template_residues`
              (which will just restrain all torsions to their current angles).
        * restrain_backbone:
            - if `True`, all phi, psi and omega dihedrals in
              `restrained_residues` hat  exist in both `restrained_residues`
              and `template_residues` will be restrained to the angles in
              `template_residues`
        * restrain_sidechains:
            - if `True`, all chi dihedrals in `restrained_residues` that  exist
              in both `restrained_residues` and `template_residues` will be
              restrained to the angles in `template_residues`
        * kappa:
            - can be thought of as (approximately) the inverse variance of the
              well within which the restraint will be felt. For example,
              kappa=10 corresponds to a standard deviation of
              :math:`\sqrt{\frac{1}{10}}` = 0.316 radians or about 18 degrees.
              A torsion restrained with this kappa will begin to "feel" the
              restraint once it comes within about two standard deviations of
              the target angle.
        * alpha:
            - sets the rate of fall-off of the potential outside the restraint
              well. For alpha==0, the potential quickly falls to a constant
              (zero force). Increasing alpha increases the force felt outside
              the core well region.
        * spring_constant:
            - strength of each restraint, in :math:`kJ mol^{-1} rad^{-2}`
        * identical_sidechains_only:
            - if True, restraints will only be applied to a sidechain if it is
              the same amino acid as the template.
        * adjust_for_confidence:
            - for use when the reference model is a predicted structure. If True, 
              the parameters of each individual restraint will be adjusted 
              according to the predicted confidence in local conformation (less
              confident restraints will be both weaker and "fuzzier")
        * confidence_type:
            - the confidence metric used in the prediction. Currently only "plddt"
              (predicted local distance difference test) is supported, and is 
              expected to be provided in the B-factor column of the reference 
              model as for models fetched from the AlphaFold2-EBI database. Scores
              should be on either a 0-100 or 0-1 scale.
    '''
    if not restrained_residues or not len(restrained_residues):
        raise UserError('No residues specified to restrain!')
    if not template_residues or not len(template_residues):
        raise UserError('No template residues specified')
    from chimerax.isolde import session_extensions as sx
    template_us = template_residues.unique_structures
    if len(template_us) != 1:
        raise UserError('Template residues must be from a single model!')
    template_model = template_us[0]
    restrained_us = restrained_residues.unique_structures
    if len(restrained_us) != 1:
        raise UserError('Restrained residues must be from a single model!')
    restrained_model = restrained_us[0]
    from ..atomic.util import correct_pseudosymmetric_sidechain_atoms
    correct_pseudosymmetric_sidechain_atoms(session, restrained_residues)
    if template_model != restrained_model:
        correct_pseudosymmetric_sidechain_atoms(session, template_residues)
    if adjust_for_confidence:
        if confidence_type=='plddt':
            confidence_multiplier=plddt_multiplier(template_model.atoms)


    tdm = sx.get_proper_dihedral_mgr(session)
    rdrm = sx.get_adaptive_dihedral_restraint_mgr(restrained_model)
    pdrm = sx.get_proper_dihedral_restraint_mgr(restrained_model)
    backbone_names = ('phi','psi')
    sidechain_names = ('chi1','chi2','chi3','chi4')
    names = []
    if restrain_backbone:
        names.extend(backbone_names)
    if restrain_sidechains:
        names.extend(sidechain_names)

    def apply_restraints(trs, rrs, adjust_for_confidence=False, confidence_type='plddt'):
        for name in names:
            for tr, rr in zip(trs, rrs):
                if 'chi' in name and identical_sidechains_only and tr.name != rr.name:
                    continue
                td = tdm.get_dihedral(tr, name)
                rdr = rdrm.add_restraint_by_residue_and_name(rr, name)
                if td and rdr:
                    if adjust_for_confidence:
                        if confidence_type=='plddt':
                            alpha_mult, k_mult, kappa_add = _adjust_torsion_restraint_terms_for_plddt(td.atoms.bfactors.min()*confidence_multiplier)
                        else:
                            raise RuntimeError(f'Unknown confidence type "{confidence_type}"!')
                        if alpha_mult == 0:
                            continue
                    else:
                        alpha_mult = 1
                        k_mult = 1
                        kappa_add = 0
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
                    rdr.spring_constant = spring_constant*k_mult
                    rdr.kappa = min(kappa+kappa_add, KAPPA_MAX)
                    rdr.alpha = alpha*alpha_mult
                    rdr.enabled = True
        if restrain_backbone:
            # For omega dihedrals we really want to stick with the standard
            # proper dihedral restraints, but we *don't* want to blindly
            # restrain them to the template. Rather, we want to set the target
            # angles to either 0 or pi.
            import numpy
            dmask = numpy.logical_and(
                tdm.residues_have_dihedral(trs, 'omega'),
                tdm.residues_have_dihedral(rrs, 'omega')
                )
            valid_trs = trs[dmask]
            valid_rrs = rrs[dmask]
            # Filter out residues that are cis-pro in the template but non-Pro
            # in the target
            tds = tdm.get_dihedrals(valid_trs, 'omega')
            from math import radians
            pro_mut_mask = numpy.logical_and(
                valid_trs.names=='PRO',
                valid_rrs.names!='PRO'
            )
            cispro_mask = numpy.logical_not(numpy.logical_and(
                pro_mut_mask, numpy.abs(tds.angles) < radians(30)
            ))
            valid_trs = valid_trs[cispro_mask]
            valid_rrs = valid_rrs[cispro_mask]

            tds = tdm.get_dihedrals(valid_trs, 'omega')
            ors = pdrm.add_restraints_by_residues_and_name(valid_rrs, 'omega')
            targets = tds.angles
            from math import radians, pi
            mask = numpy.abs(targets) > radians(30)
            ors.targets = numpy.ones(len(mask))*pi*mask
            ors.enableds = True
            ors.displays = False

    if template_residues==restrained_residues:
        apply_restraints(template_residues, restrained_residues, adjust_for_confidence, confidence_type)

    else:
        # Just do it based on sequence alignment for now.
        # TODO: work out an improved weighting scheme based on similarity of
        # environment for each residue.
        trs, rrs = sequence_align_all_residues(
            session, [template_residues], [restrained_residues]
        )
        apply_restraints(trs, rrs, adjust_for_confidence, confidence_type)

def adjust_torsion_restraint_terms_for_confidence(score, confidence_type):
    if confidence_type=='plddt':
        return _adjust_torsion_restraint_terms_for_plddt(score)

def _adjust_torsion_restraint_terms_for_plddt(score):
    if score < 0.5:
        return 0,0,0
    from math import exp
    alpha_multiplier=exp(-10*(score-1)**2)
    spring_constant_multiplier= exp(-6*(score-1)**2)
    kappa_add = 2400*(exp(-9.5*score)-exp(-9.5))
    return alpha_multiplier, spring_constant_multiplier, kappa_add


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
    from chimerax.geometry import find_close_points, distance
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


def restrain_small_ligands(model, distance_cutoff=4, heavy_atom_limit=3, spring_constant=5000,
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

    from chimerax.geometry import find_close_points, distance
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
            num_drs = 0
            for fi in found_i:
                if fi in r_indices:
                    continue
                dr = drm.add_restraint(ra, all_heavy_atoms[fi])
                dr.spring_constant = spring_constant
                dr.target=distance(rc, all_heavy_coords[fi])
                dr.enabled = True
                num_drs += 1
                # applied_drs = True
        if num_drs < 3:
            # Really need at least 3 distance restraints (probably 4, actually,
            # but we don't want to be *too* restrictive) to be stable in 3D
            # space. If we have fewer than that, add position restraints to be
            # sure.
            prs = prm.add_restraints(r_heavy_atoms)
            prs.targets = prs.atoms.coords
            prs.spring_constants = spring_constant
            prs.enableds = True


def sequence_align_all_residues(session, residues_a, residues_b):
    '''
    Peform a sequence alignment on each pair of chains, and concatenate the
    aligned residues into a single "super-alignment" (a pair of Residues
    objects with 1:1 correspondence between residues at each index). Each
    :class:`Residues` object in the arguments should contain residues from a
    single chain, and the chains in :param:`residues_a` and `residues_b` should
    be matched.

    Arguments:

        * session:
            - the ChimeraX session object
        * residues_a:
            - a list of ChimeraX :class:`Residues` objects
        * residues_b:
            - a list of ChimeraX :class:`Residues` objects
    '''
    from chimerax.match_maker.match import defaults, align, check_domain_matching
    alignments = ([],[])
    dssp_cache={}
    for ra, rb in zip(residues_a, residues_b):
        uca = ra.unique_chains
        ucb = rb.unique_chains
        if not (len(uca)==1 and len(ucb)==1):
            raise TypeError('Each residue selection must be from a single chain!')
        match, ref = [check_domain_matching([ch],dr)[0] for ch, dr in
            ((uca[0], ra), (ucb[0], rb))]
        score, s1, s2 = align(
            session, match, ref, defaults['matrix'],
            'nw', defaults['gap_open'],
            defaults['gap_extend'], dssp_cache
        )
        for i, (rr, mr) in enumerate(zip(s1, s2)):
            if mr=='.' or rr=='.':
                continue
            ref_res = s1.residues[s1.gapped_to_ungapped(i)]
            match_res = s2.residues[s2.gapped_to_ungapped(i)]
            if not ref_res or not match_res:
                continue
            alignments[0].append(ref_res)
            alignments[1].append(match_res)
    from chimerax.atomic import Residues
    return [Residues(a) for a in alignments]

def paired_principal_atoms(aligned_residues):
    from chimerax.atomic import Atoms
    result = [Atoms(aa) for aa in zip(*((a1, a2) for a1, a2 in zip(
        aligned_residues[0].principal_atoms, aligned_residues[1].principal_atoms
    ) if a1 is not None and a2 is not None
    ))]
    if not len(result):
        return [Atoms(), Atoms()]
    return result




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

def restrain_atom_distances_to_template(session, template_residues, restrained_residues,
    protein=True, nucleic=True, custom_atom_names=[],
    distance_cutoff=8, alignment_cutoff=5, well_half_width = 0.1,
    kappa = 10, tolerance = 0.025, fall_off = 2, display_threshold=None,
    adjust_for_confidence=False, use_coordinate_alignment=True, confidence_type='pae', pae_matrix=None):
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
              value after alignment will not be restrained. Ignored if `use_coordinate_alignment`
              is `False`.
        * well_half_width (default = 0.1):
            - distance range (as a fraction of the square root of target distance) within which
              the restraint will behave like a normal harmonic restraint.
              The applied force will gradually taper off for any restraint
              deviating from (target + tolerance) by more than this amount.
        * kappa (default = 10):
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
        * fall_off (default = 2):
            - Sets the rate at which the energy function will fall off when the
              distance deviates strongly from the target, as a function of the
              target distance. The exponent on the energy term at large
              deviations from the target distance will be set as
              :math:`\alpha = -1 -\text{fall\_off} ln(\text{target})`. In other
              words, long-distance restraints are treated as less confident than
              short-distance ones.
        * display_threshold (default = 0):
            - deviation from (target +- tolerance) as a fraction of
              :attr:`well_half_width` below which restraints will be hidden.
        * use_coordinate_alignment (default = True):
            - if True, reference and template residues will be matched by progressively breaking down
              the models into groups that approximately align as rigid bodies. This is usually the 
              preferable approach, but fails in cases where the working model is badly wrong (e.g. has 
              large register errors). If False, residues will be matched strictly by sequence alignment.
              This may be preferable when restraining against an AlphaFold model.
        * adjust_for_confidence (default = False):
            - if true, interpret B-factors of the template atoms as a confidence score, the 
              definition of which is controlled by :attr:`confidence_type`, and adjust 
              :attr:`kappa`, :attr:`tolerance` and :attr:`fall_off` according to the mean 
              confidence for each restrained atom pair.
        * confidence_type (default = "pae"):
            - the type of confidence score used. Current options are "pae" (predicted aligned error) and
              "plddt" (predicted local distance difference test). For the "pae" option, multi-chain template 
              models are NOT currently supported.
        * pae_matrix (default = None):
            - used if adjust_for_confidence is True and confidence_type is "pae". If the reference model
              was downloaded from the AlphaFold database, leave this argument as None and the relevant PAE
              matrix will be automatically fetched. Otherwise, the matrix should be a 2D Numpy array with entry (i,j) 
              equal to the PAE of residue number i relative to residue number j. 

    '''
    from chimerax.std_commands.align import IterationError
    from chimerax.isolde import session_extensions as sx
    if not protein and not nucleic and not len(custom_atom_names):
        raise UserError('Nothing to restrain!')
    # if len(template_residues) != len(restrained_residues):
    #     raise TypeError('Template and restrained residue arrays must be the same length!')
    for rrs in restrained_residues:
        if len(rrs) == 0:
            raise UserError('No residues specified to restrain!')
        restrained_us = rrs.unique_structures
        if len(restrained_us) != 1:
            raise UserError('Restrained residues must be from a single model!')
    for trs in template_residues:
        if len(trs) == 0:
            raise UserError('No template residues specified!')
        template_us = trs.unique_structures
        if len(template_us) != 1:
            raise UserError('Template residues must be from a single model! '
                'Residues are {} in {}'.format(trs.numbers, ','.join(s.id_string for s in trs.structures)))
    restrained_model = restrained_us[0]
    from ..atomic.util import correct_pseudosymmetric_sidechain_atoms
    correct_pseudosymmetric_sidechain_atoms(session, restrained_model.residues)
    for tu in template_us:
        if tu != restrained_model:
            correct_pseudosymmetric_sidechain_atoms(session, tu.residues)
    log = restrained_model.session.logger
    if adjust_for_confidence and confidence_type not in ('plddt','pae'):
        raise UserError('confidence_type must be one of ("plddt", "pae")!')

    if adjust_for_confidence:
        if confidence_type=='pae':
            chains = set()
            for trs in template_residues:
                chains.update(set(trs.unique_chain_ids))
            if len(chains) != 1:
                raise UserError('Weighting according to PAE is currently only supported for single-chain templates!')
            if pae_matrix is None:
                from chimerax.isolde.reference_model.alphafold import alphafold_id, fetch_alphafold_pae
                tm = template_residues[0].unique_structures[0]
                aid = alphafold_id(tm)
                if aid is None:
                    raise UserError(f'Template model #{tm.id_string} does not appear to be from the AlphaFold server. You will need to provide a PAE matrix separately.')
                uniprot_id = aid.split('-')[1]
                try:
                    pae_matrix = fetch_alphafold_pae(session, uniprot_id)
                except:
                    raise UserError(f'Failed to fetch the PAE matrix for template model #{tm.id_string}!')
        elif confidence_type=='plddt':
            confidence_multiplier = plddt_multiplier(template_us.atoms)



    adrm = sx.get_adaptive_distance_restraint_mgr(restrained_model)
    if display_threshold is not None:
        adrm.display_threshold = display_threshold
    from chimerax.geometry import find_close_points, distance
    from chimerax.atomic import concatenate

    atom_names = []
    if protein:
        atom_names.extend(['CA', 'CB', 'CG', 'CG1', 'OG', 'OG1'])
    if nucleic:
        atom_names.extend(["OP1", "OP2", "C4'", "C2'", "O2", "O4", "N4", "N2", "O6", "N1", "N6", "N9"])
    atom_names.extend(custom_atom_names)

    import numpy
    atom_names = set(atom_names)

    def apply_restraints(trs, rrs, adjust_for_confidence, confidence_type):
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
        from math import sqrt
        for i, ra1 in enumerate(restrained_as):
            query_coord = numpy.array([template_coords[i]])
            indices = find_close_points(query_coord, template_coords, distance_cutoff)[1]
            indices = indices[indices !=i]
            for ind in indices:
                ra2 = restrained_as[ind]
                if ra1.residue == ra2.residue:
                    continue
                if adjust_for_confidence:
                    if confidence_type=='plddt':
                        scores = [template_as[i].bfactor*confidence_multiplier, template_as[ind].bfactor*confidence_multiplier]
                    elif confidence_type=='pae':
                        scores = [pae_matrix[template_as[i].residue.number-1,template_as[ind].residue.number-1]]
                    kappa_adj, tol_adj, falloff_adj = adjust_distance_restraint_terms_by_confidence(scores, confidence_type)
                    if kappa_adj == 0:
                        continue
                else:
                    kappa_adj = tol_adj = 1
                    falloff_adj = 0
                try:
                    dr = adrm.add_restraint(ra1, ra2)
                except ValueError:
                    continue
                dist = distance(query_coord[0], template_coords[ind])
                dr.tolerance = tolerance * dist * tol_adj
                dr.target = dist
                dr.c = max(sqrt(dist)*well_half_width, 0.1)
                #dr.effective_spring_constant = spring_constant
                dr.kappa = kappa * kappa_adj
                from math import log
                dr.alpha = -1 - fall_off * log((max(dist-1,1))) - falloff_adj
                dr.enabled = True

    if all(trs == rrs for trs, rrs in zip(template_residues, restrained_residues)):
        # If the template is identical to the model, we can just go ahead and restrain

        return [apply_restraints(trs, rrs, adjust_for_confidence, confidence_type) for trs, rrs in zip(template_residues, restrained_residues)]
    else:
        if not use_coordinate_alignment:
            atrs, arrs = sequence_align_all_residues(session, template_residues, restrained_residues)
            return apply_restraints(atrs, arrs, adjust_for_confidence, confidence_type)
        else:
            # If we're not simply restraining the model to itself, we need to do an
            # alignment to get a matching pair of residue sequences
            tpa, rpa = paired_principal_atoms(sequence_align_all_residues(
                session, template_residues, restrained_residues
            ))
            from chimerax.std_commands import align
            while len(tpa) >3 and len(rpa) >3:
                try:
                    ta, ra, *_ = align.align(session, tpa, rpa, move=False,
                        cutoff_distance=alignment_cutoff)
                except IterationError:
                    log.info(('No further alignments of 3 or more residues found.'))
                    break
                tr, rr = [ta.residues, ra.residues]
                apply_restraints(tr, rr, adjust_for_confidence, confidence_type)
                tpa = tpa.subtract(ta)
                rpa = rpa.subtract(ra)



        # while len(template_residues) != 0 and len(restrained_residues) != 0:
        #     found_trs = []
        #     found_rrs = []
        #     remain_trs = []
        #     remain_rrs = []
        #     for trs, rrs in zip(template_residues, restrained_residues):
        #         try:
        #             ftrs, frrs = _get_template_alignment(
        #                 trs, rrs, cutoff_distance = alignment_cutoff,
        #                 overlay_template = False, always_raise_errors = True
        #             )
        #         except IterationError:
        #             log.info(("No sufficient alignment found for {} residues "
        #                 "in template chain {} and {} residues in restrained "
        #                 "model chain {}").format(len(trs), trs.chain_ids[0],
        #                 len(rrs), rrs.chain_ids[0]))
        #             continue
        #         found_trs.append(ftrs)
        #         found_rrs.append(frrs)
        #         rtrs = trs.subtract(ftrs)
        #         rrrs = rrs.subtract(frrs)
        #         if len(rtrs) > 3 and len(rrrs) > 3:
        #             remain_trs.append(rtrs)
        #             remain_rrs.append(rrrs)
        #     template_residues = remain_trs
        #     restrained_residues = remain_rrs
        #
        #     if len(found_trs):
        #         trs = concatenate(found_trs)
        #         rrs = concatenate(found_rrs)
        #         apply_restraints(trs, rrs)

def plddt_multiplier(atoms):
    '''
    Different implementations of the AlphaFold2 algorithm provide pLDDT in the B-factor column on 
    either a 0-1 or 0-100 range. For normalisation purposes, we assume that if any atom has a 
    B-factor over 1 then the range is 0-100. 
    '''
    if atoms.bfactors.max() <= 1:
        return 1.0
    return 0.01

def adjust_distance_restraint_terms_by_confidence(confidence_scores, confidence_type='plddt'):
    if confidence_type=='plddt':
        return adjust_distance_restraint_terms_by_predicted_lddt(confidence_scores)
    elif confidence_type=='pae':
        return adjust_distance_restraint_terms_by_pae(confidence_scores[0])
    from chimerax.core.errors import UserError
    raise UserError(f'No confidence type "{confidence_type}"!')


def adjust_distance_restraint_terms_by_predicted_lddt(confidence_scores):
    from math import exp
    # Use the least-confident position for weighting purposes
    c = min(confidence_scores)
    if c<0.5:
        return 0,0,0
    kappa_adj = exp(-5*(c-1)**2)
    tol_adj = exp(2*(c-1)**2)
    falloff_adj = exp(-3*(c-1))-1
    return kappa_adj, tol_adj, falloff_adj

def adjust_distance_restraint_terms_by_pae(pae):
    if pae > 4:
        return 0,0,0
    from math import exp
    kappa_adj = exp(-0.1*pae**2)
    tol_adj = exp((pae/4)**2)
    falloff_adj = (pae-0.2)/2
    return kappa_adj, tol_adj, falloff_adj





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
