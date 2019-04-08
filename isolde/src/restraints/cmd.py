# @Author: Tristan Croll <tic20>
# @Date:   05-Apr-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 05-Apr-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll

def restrain_distances(session, atoms, template_atoms=None, **kwargs):
    valid_args = set((
        'protein', 'nucleic', 'custom_atom_names', 'distance_cutoff',
        'alignment_cutoff', 'well_half_width', 'kappa', 'tolerance', 'fall_off'
    ))
    log = session.logger

    from . import restraint_utils

    model_residues = atoms.unique_residues
    if template_atoms is None:
        template_residues = model_residues
    else:
        template_residues = template_atoms.unique_residues
    args = {kw: arg for kw, arg in kwargs.items() if kw in valid_args and arg is not None}
    restraint_utils.restrain_atom_distances_to_template(template_residues,
        model_residues, **kwargs)

def release_adaptive_distance_restraints(session, atoms,
        longer_than=None, strained_only=False, stretch_limit=1.5, compression_limit=0.5):
    from chimerax.isolde import session_extensions as sx
    for m in atoms.unique_structures:
        adrm = sx.get_adaptive_distance_restraint_mgr(m)
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
        intra_only=False, kappa=None, well_half_width=None,
        tolerance=None, fall_off=None):
    log = session.logger
    from chimerax.isolde import session_extensions as sx
    for m in atoms.unique_structures:
        adrm = sx.get_adaptive_distance_restraint_mgr(m)
        if intra_only:
            adrs = adrm.intra_restraints(atoms)
        else:
            adrs = adrm.atoms_restraints(atoms)
        if kappa:
            adrs.kappas = kappa
        targets = adrs.targets
        if well_half_width:
            adrs.cs = well_half_width * targets
        if tolerance:
            adrs.tolerances = tolerance*targets
        if fall_off:
            from math import log
            import numpy
            alphas = numpy.empty(targets.shape, numpy.double)
            dmask = (targets < 1)
            alphas[mask] = -2
            dmask = numpy.logical_not(dmask)
            alphas[dmask] = -2-fall_off*log(targets[dmask])
            adrs.alphas = alphas



def register_isolde_restrain(logger):
    from chimerax.core.commands import (
        register, CmdDesc,
        AtomSpecArg, ModelArg,
        FloatArg, IntArg, BoolArg, StringArg, NoArg,
        ListOf, EnumOf, RepeatOf
    )
    from chimerax.atomic import AtomsArg

    def register_isolde_restrain_distances():
        desc = CmdDesc(
            synopsis = 'Restrain interatomic distances to current values or a template.',
            required = [
                ('atoms', AtomsArg),
            ],
            keyword = [
                ('template_atoms', AtomsArg),
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
    def register_isolde_release_distances():
        desc = CmdDesc(
            synopsis = 'Release overly-strained adaptive distance restraints',
            required = [('atoms', AtomsArg)],
            keyword = [
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
                ('intra_only', BoolArg),
                ('kappa', FloatArg),
                ('well_half_width', FloatArg),
                ('tolerance', FloatArg),
                ('fall_off', FloatArg)
            ]
        )
        register('isolde adjust distances', desc, adjust_adaptive_distance_restraints, logger=logger)
    register_isolde_restrain_distances()
    register_isolde_release_distances()
    register_isolde_adjust_distances()
