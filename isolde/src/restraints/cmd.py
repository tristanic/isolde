# @Author: Tristan Croll <tic20>
# @Date:   05-Apr-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 19-Nov-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll

from chimerax.core.errors import UserError

def restrain_distances(session, atoms, template_atoms=None, per_chain=False, **kwargs):
    valid_args = set((
        'protein', 'nucleic', 'custom_atom_names', 'distance_cutoff',
        'alignment_cutoff', 'well_half_width', 'kappa', 'tolerance', 'fall_off',
        'display_threshold'
    ))
    log = session.logger

    from . import restraint_utils

    if template_atoms is not None and len(atoms) != len(template_atoms):
        log.warning('You must provide one template selection for each restrained selection!')
        return

    if per_chain and (template_atoms is not None):
        log.warning('perChain argument is only valid when restraining atoms to '
            'their initial conformation, not when restraining to a template.')
        return

    model_residues = [mas.unique_residues for mas in atoms]
    if per_chain:
        split_residues=[]
        for residues in model_residues:
            chain_ids = residues.unique_chain_ids
            for cid in chain_ids:
                split_residues.append(residues[residues.chain_ids==cid])
        model_residues = split_residues
    if template_atoms is None:
        template_residues = model_residues
    else:
        template_residues = [tas.unique_residues for tas in template_atoms]
        for mr, tr in zip(model_residues, template_residues):
            if len(mr.unique_chains) != 1 or len(tr.unique_chains) != 1:
                raise UserError('Each atom selection should be from a single chain!')

    args = {kw: arg for kw, arg in kwargs.items() if kw in valid_args and arg is not None}
    restraint_utils.restrain_atom_distances_to_template(session, template_residues,
        model_residues, **kwargs)

def release_adaptive_distance_restraints(session, atoms,
        internal_only=False, external_only=False,
        longer_than=None, strained_only=False,
        stretch_limit=1.2, compression_limit=0.8):
    log = session.logger
    if internal_only and external_only:
        log.warning('Cannot specify both internal only and external only!')
        return
    from chimerax.isolde import session_extensions as sx
    for m in atoms.unique_structures:
        adrm = sx.get_adaptive_distance_restraint_mgr(m)
        if internal_only:
            adrs = adrm.intra_restraints(atoms)
        elif external_only:
            adrs = adrm.atoms_restraints(atoms).subtract(adrm.intra_restraints(atoms))
        else:
            adrs = adrm.atoms_restraints(atoms)
        if longer_than:
            adrs = adrs[adrs.distances > longer_than]
        if not strained_only:
            adrs.enableds = False
        else:
            if stretch_limit is not None:
                adrs[adrs.distances/(adrs.targets+adrs.tolerances)>stretch_limit].enableds = False
            if compression_limit is not None:
                non_zeros = adrs[adrs.targets-adrs.tolerances > 0]
                non_zeros[non_zeros.distances/(non_zeros.targets-non_zeros.tolerances)<compression_limit].enableds = False

def adjust_adaptive_distance_restraints(session, atoms,
        internal_only=False, external_only=False, kappa=None, well_half_width=None,
        tolerance=None, fall_off=None, display_threshold=None):
    log = session.logger
    if internal_only and external_only:
        log.warning('Cannot specify both internal only and external only!')
        return
    from chimerax.isolde import session_extensions as sx
    for m in atoms.unique_structures:
        adrm = sx.get_adaptive_distance_restraint_mgr(m)
        if internal_only:
            adrs = adrm.intra_restraints(atoms)
        elif external_only:
            adrs = adrm.atoms_restraints(atoms).subtract(adrm.intra_restraints(atoms))
        else:
            adrs = adrm.atoms_restraints(atoms)
        if kappa is not None:
            adrs.kappas = kappa
        targets = adrs.targets
        if well_half_width is not None:
            adrs.cs = well_half_width * targets
        if tolerance is not None:
            adrs.tolerances = tolerance*targets
        if fall_off is not None:
            import numpy
            alphas = numpy.empty(targets.shape, numpy.double)
            dmask = (targets < 1)
            alphas[dmask] = -2
            dmask = numpy.logical_not(dmask)
            alphas[dmask] = -2-fall_off*numpy.log(targets[dmask])
            adrs.alphas = alphas
        if display_threshold is not None:
            adrm.display_threshold = display_threshold

def restrain_single_distance(session, atoms, min_dist, max_dist,
        strength=20, well_half_width=None, confidence=-2):
    if len(atoms)!=2:
        raise UserError('Selection must specify exactly two atoms!')
    atom1, atom2 = atoms
    if not atom1.structure == atom2.structure:
        raise UserError('Both atoms must be from the same model!')
    target = (min_dist + max_dist) / 2
    tolerance = (max_dist - min_dist) / 2
    if well_half_width is None:
        well_half_width = min(target/5, 2.0)
    from .restraint_utils import restrain_atom_pair_adaptive_distance
    restrain_atom_pair_adaptive_distance(atom1, atom2, target, tolerance,
        strength, well_half_width, confidence)

def restrain_ligands(session, models, distance_cutoff=4, max_heavy_atoms=3,
        spring_constant=5000, bond_to_carbon=False):
    if not models:
        raise UserError('Must specify at least one atomic structure!')
    from chimerax.isolde.restraints.restraint_utils import restrain_small_ligands
    for m in models:
        restrain_small_ligands(m, distance_cutoff=distance_cutoff,
            heavy_atom_limit=max_heavy_atoms, spring_constant=spring_constant,
            bond_to_carbon=bond_to_carbon)

def _dihedral_names(backbone, sidechains):
    backbone_names = ('phi', 'psi')
    sidechain_names = ('chi1', 'chi2', 'chi3', 'chi4')
    names = []
    if backbone:
        names.extend(backbone_names)
    if sidechains:
        names.extend(sidechain_names)
    return names

def restrain_torsions(session, residues, template_residues=None,
        backbone=True, sidechains=True, angle_range=None, alpha=None,
        spring_constant=250, identical_sidechains_only=True, adjust_for_confidence=False,
        confidence_type='plddt'):
    if angle_range is None:
        if adjust_for_confidence:
            angle_range = 150
        else:
            angle_range = 60
    if alpha is None:
        if adjust_for_confidence:
            alpha = 0.5
        else:
            alpha = 0.2

    logger = session.logger
    if not residues:
        raise UserError('Must specify residues to restrain')
    from chimerax.atomic import Residue
    pres = residues[residues.polymer_types==Residue.PT_AMINO]
    if not len(pres):
        raise UserError('This command is only applicable to amino acids, but the '
            'selection contains none.')
    if len(pres) != len(residues):
        logger.warning('The "isolde restrain torsions" command only applies to '
            'protein chains. Other residues have been ignored.')
    if template_residues:
        tres = template_residues[template_residues.polymer_types==Residue.PT_AMINO]
    else:
        tres = pres

    from math import radians
    from ..molobject import AdaptiveDihedralRestraintMgr
    kappa = AdaptiveDihedralRestraintMgr.angle_range_to_kappa(radians(angle_range))

    from .restraint_utils import restrain_torsions_to_template
    restrain_torsions_to_template(session, tres, pres,
        restrain_backbone=backbone,
        restrain_sidechains=sidechains, kappa=kappa, alpha=alpha,
        spring_constant=spring_constant,
        identical_sidechains_only=identical_sidechains_only,
        adjust_for_confidence=adjust_for_confidence, confidence_type=confidence_type)

def adjust_torsions(session, residues, backbone=True,
        sidechains=True, angle_range=None, alpha=None, spring_constant=None):
    from chimerax.isolde.session_extensions import get_adaptive_dihedral_restraint_mgr
    if angle_range is not None:
        from ..molobject import AdaptiveDihedralRestraintMgr
        from math import radians
        kappa = AdaptiveDihedralRestraintMgr.angle_range_to_kappa(radians(angle_range))
    names = _dihedral_names(backbone, sidechains)
    for m in residues.unique_structures:
        apdrm = get_adaptive_dihedral_restraint_mgr(m, create=False)
        if apdrm is None:
            continue
        for name in names:
            apdrs = apdrm.get_restraints_by_residues_and_name(residues, name)
            if angle_range is not None:
                apdrs.kappas = kappa
            if spring_constant is not None:
                apdrs.spring_constants = spring_constant
            if alpha is not None:
                apdrs.alphas = alpha

def release_torsions(session, residues, backbone=True, sidechains=True):
    from chimerax.isolde.session_extensions import get_adaptive_dihedral_restraint_mgr
    names = _dihedral_names(backbone, sidechains)
    for m in residues.unique_structures:
        apdrm = get_adaptive_dihedral_restraint_mgr(m, create=False)
        if apdrm is None:
            continue
        for name in names:
            apdrm.get_restraints_by_residues_and_name(residues, name).enableds=False



def register_isolde_restrain(logger):
    from chimerax.core.commands import (
        register, CmdDesc,
        AtomSpecArg, ModelArg,
        FloatArg, IntArg, BoolArg, StringArg, NoArg,
        ListOf, EnumOf, RepeatOf
    )
    from chimerax.atomic import AtomsArg, AtomicStructuresArg, ResiduesArg

    def register_isolde_restrain_distances():
        desc = CmdDesc(
            synopsis = 'Restrain interatomic distances to current values or a template.',
            required = [
                ('atoms', ListOf(AtomsArg)),
            ],
            keyword = [
                ('template_atoms', ListOf(AtomsArg)),
                ('protein', BoolArg),
                ('nucleic', BoolArg),
                ('custom_atom_names', ListOf(StringArg)),
                ('per_chain', BoolArg),
                ('distance_cutoff', FloatArg),
                ('use_coordinate_alignment', BoolArg),
                ('alignment_cutoff', FloatArg),
                ('well_half_width', FloatArg),
                ('kappa', FloatArg),
                ('tolerance', FloatArg),
                ('fall_off', FloatArg),
                ('display_threshold', FloatArg),
                ('adjust_for_confidence', BoolArg),
                ('confidence_type', EnumOf(['plddt','pae']))
            ]
        )
        register('isolde restrain distances', desc, restrain_distances, logger=logger)
    def register_isolde_restrain_single_distance():
        desc = CmdDesc(
            synopsis = 'Create a single adaptive distance restraint',
            required = [
                ('atoms', AtomsArg),
                ('min_dist', FloatArg),
                ('max_dist', FloatArg),
            ],
            keyword = [
                ('strength', FloatArg),
                ('well_half_width', FloatArg),
                ('confidence', FloatArg),
            ]
        )
        register('isolde restrain single distance', desc, restrain_single_distance, logger=logger)
    def register_isolde_release_distances():
        desc = CmdDesc(
            synopsis = 'Release selected adaptive distance restraints',
            required = [('atoms', AtomsArg)],
            keyword = [
                ('internal_only', BoolArg),
                ('external_only', BoolArg),
                ('longer_than', FloatArg),
                ('strained_only', BoolArg),
                ('stretch_limit', FloatArg),
                ('compression_limit', FloatArg)
            ]
        )
        register('isolde release distances', desc, release_adaptive_distance_restraints, logger=logger)
    def register_isolde_adjust_distances():
        desc = CmdDesc(
            synopsis = 'Adjust the behaviour of existing adaptive distance restraints',
            required = [('atoms', AtomsArg)],
            keyword = [
                ('internal_only', BoolArg),
                ('external_only', BoolArg),
                ('kappa', FloatArg),
                ('well_half_width', FloatArg),
                ('tolerance', FloatArg),
                ('fall_off', FloatArg),
                ('display_threshold', FloatArg)
            ]
        )
        register('isolde adjust distances', desc, adjust_adaptive_distance_restraints, logger=logger)
    def register_isolde_restrain_ligands():
        desc = CmdDesc(
            synopsis = ('Restrain small ligands to their surroundings and/or '
                'current positions. If less than four distance restraints can '
                'be made for a given residue, they will be supplemented with '
                'position restraints.'),
            required = [
                ('models', AtomicStructuresArg),
            ],
            keyword = [
                ('distance_cutoff', FloatArg),
                ('max_heavy_atoms', IntArg),
                ('spring_constant', FloatArg),
                ('bond_to_carbon', BoolArg),
            ]
        )
        register('isolde restrain ligands', desc, restrain_ligands, logger=logger)

    def register_isolde_restrain_torsions():
        desc = CmdDesc(
            synopsis = ('Restrain amino acid backbone and/or sidechain torsions '
                'to their current values or to those of a template model of '
                'identical or homologous sequence.'),
            required = [
                ('residues', ResiduesArg),
            ],
            keyword = [
                ('template_residues', ResiduesArg),
                ('backbone', BoolArg),
                ('sidechains', BoolArg),
                ('angle_range', FloatArg),
                ('alpha', FloatArg),
                ('spring_constant', FloatArg),
                ('identical_sidechains_only', BoolArg),
                ('adjust_for_confidence', BoolArg),
            ]
        )
        register('isolde restrain torsions', desc, restrain_torsions, logger=logger)

    def register_isolde_adjust_torsions():
        desc = CmdDesc(
            synopsis = ('Adjust the strength and/or angle ranges of torsion '
                'restraints previously applied using "isolde restrain torsions".'
                ),
            required = [
                ('residues', ResiduesArg),
            ],
            keyword = [
                ('backbone', BoolArg),
                ('sidechains', BoolArg),
                ('angle_range', FloatArg),
                ('alpha', FloatArg),
                ('spring_constant', FloatArg),
            ]
        )
        register('isolde adjust torsions', desc, adjust_torsions, logger=logger)

    def register_isolde_release_torsions():
        desc = CmdDesc(
            synopsis=('Release a set of torsion restraints previously applied '
                'using "isolde restrain torsions".'),
            required = [
                ('residues', ResiduesArg),
            ],
            keyword = [
                ('backbone', BoolArg),
                ('sidechains', BoolArg),
            ]
        )
        register('isolde release torsions', desc, release_torsions, logger=logger)


    register_isolde_restrain_distances()
    register_isolde_restrain_single_distance()
    register_isolde_release_distances()
    register_isolde_adjust_distances()
    register_isolde_restrain_ligands()
    register_isolde_restrain_torsions()
    register_isolde_adjust_torsions()
    register_isolde_release_torsions()
