# @Author: Tristan Croll <tic20>
# @Date:   05-Apr-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 12-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll

from chimerax.core.errors import UserError

def restrain_distances(session, atoms, template_atoms=None, **kwargs):
    valid_args = set((
        'protein', 'nucleic', 'custom_atom_names', 'distance_cutoff',
        'alignment_cutoff', 'well_half_width', 'kappa', 'tolerance', 'fall_off'
    ))
    log = session.logger

    from . import restraint_utils

    if template_atoms is not None and len(atoms) != len(template_atoms):
        log.warning('You must provide one template selection for each restrained selection!')
        return

    model_residues = [mas.unique_residues for mas in atoms]
    if template_atoms is None:
        template_residues = model_residues
    else:
        template_residues = [tas.unique_residues for tas in template_atoms]
        for mr, tr in zip(model_residues, template_residues):
            if len(mr.unique_chains) != 1 or len(tr.unique_chains) != 1:
                raise UserError('Each atom selection should be from a single chain!')

    args = {kw: arg for kw, arg in kwargs.items() if kw in valid_args and arg is not None}
    restraint_utils.restrain_atom_distances_to_template(template_residues,
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




def register_isolde_restrain(logger):
    from chimerax.core.commands import (
        register, CmdDesc,
        AtomSpecArg, ModelArg,
        FloatArg, IntArg, BoolArg, StringArg, NoArg,
        ListOf, EnumOf, RepeatOf
    )
    from chimerax.atomic import AtomsArg, AtomicStructuresArg

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
                ('distance_cutoff', FloatArg),
                ('alignment_cutoff', FloatArg),
                ('well_half_width', FloatArg),
                ('kappa', FloatArg),
                ('tolerance', FloatArg),
                ('fall_off', FloatArg)
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


    register_isolde_restrain_distances()
    register_isolde_restrain_single_distance()
    register_isolde_release_distances()
    register_isolde_adjust_distances()
    register_isolde_restrain_ligands()
