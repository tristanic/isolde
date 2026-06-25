from chimerax.core.errors import UserError

def _check_for_isolde(session):
    if not hasattr(session, 'isolde'):
        raise UserError('ISOLDE must be started first!')
    return session.isolde

def _eligible_residues(isolde, atoms):
    m = isolde.selected_model
    if m is None:
        raise UserError('No model currently selected for ISOLDE!')
    from chimerax.atomic import Residue
    residues = atoms.unique_residues
    residues = residues[residues.polymer_types==Residue.PT_AMINO]
    residues = m.residues.intersect(residues)
    if isolde.simulation_running:
        sc = isolde.sim_manager.sim_construct
        residues = sc.mobile_residues.intersect(residues)
    return residues

def pep_flip(session, atoms):
    '''
    Attempt to flip the peptide bond N-terminal to each eligible residue in the
    selection by applying temporary restraints on their phi and psi dihedrals
    180 degrees away from the current angles. ISOLDE must be initialised, and
    atoms outside of the current working model (or outside of the mobile
    selection if a simulation *is* running) will be ignored.
    '''
    isolde = _check_for_isolde(session)
    residues = _eligible_residues(isolde, atoms)
    if not len(residues):
        return
    from .peptide_flip import Peptide_Bond_Flipper
    if not isolde.simulation_running:
        from ..cmd import isolde_sim
        isolde_sim(session, 'start', residues.atoms)
    for r in residues:
        try:
            Peptide_Bond_Flipper(isolde, r)
        except TypeError:
            continue

def cis_flip(session, atoms):
    '''
    Flip the omega (peptide plane) dihedral N-terminal to each eligible residue
    from cis to trans or vice-versa. ISOLDE must be initialised, and
    atoms outside of the current working model (or outside of the mobile
    selection if a simulation *is* running) will be ignored.
    '''
    isolde = _check_for_isolde(session)
    residues = _eligible_residues(isolde, atoms)
    if not len(residues):
        return
    if not isolde.simulation_running:
        from ..cmd import isolde_sim
        isolde_sim(session, 'start', residues.atoms)
    for r in residues:
        try:
            isolde.flip_peptide_omega(r)
        except TypeError:
            continue



def chiral_flip(session, atoms, force=False):
    '''
    Repair inverted chiral centres in the selection by rigidly rotating the
    smallest substituent 180 degrees about the centre, starting a suitable
    simulation if none is running so the (deliberately strained) result relaxes
    into the correct chiral well. ISOLDE must be initialised.

    By default only centres whose handedness is actually *wrong* are flipped: a
    centre that is already correct (even if strained) is left untouched, since
    flipping it would invert a correct centre. Centres outside the working model
    (or outside the mobile selection if a simulation is running) are ignored, and
    centres with no small rotatable substituent (e.g. fused-ring centres) are
    skipped and reported -- those require delete-and-rebuild instead.

    **force (experts only):** if ``force`` is True, the reference handedness of
    *every* selected centre is flipped to its mirror image (each restraint target
    is negated), and atoms are moved only where needed to suit: a centre that
    currently *matches* its pre-flip target is geometry-inverted (it would
    otherwise fight the flipped target), while one that is *already inverted*
    relative to its pre-flip target is left untouched (it already matches the
    flipped target). The result is the geometry most likely to energy-minimise to
    satisfy all the new targets, moving the fewest atoms. This is meant only for
    actively building novel ligands or exploring alternative stereoisomers
    ("what-if"). The target change is **session-transient and per-centre**: it
    mutates only the selected :class:`ChiralCenter` instances, not the
    CCD/``chirals.json`` definition, so reloading the model (or rebuilding the
    chiral definitions) restores the original reference stereochemistry. Use with
    care -- a force-flipped centre will read as *correct* in validation, so an
    unintended force flip can mask a genuine stereochemical error.

    **Making a force flip permanent.** Because the chiral definition is keyed on
    the residue *name*, a future session will flip the centre straight back to the
    CCD-defined hand. To prevent that, give the residue a new (unique) name -- the
    renamed residue no longer matches the CCD entry, so nothing overrides your
    chosen stereochemistry. The trade-off: a renamed residue has *no* chiral
    definition of its own, so its chirality will be **unrestrained** in future
    sessions (nothing holds it in either hand) unless you also register it as a
    custom template. Full template-registration tooling is still in progress, but
    getting close.
    '''
    import numpy
    isolde = _check_for_isolde(session)
    m = isolde.selected_model
    if m is None:
        raise UserError('No model currently selected for ISOLDE!')
    from chimerax.isolde.molobject import get_chiral_mgr
    cm = get_chiral_mgr(session)
    chirals = cm.get_chirals(atoms, create=True)
    # Restrict to centres in the working model, or the mobile selection mid-sim.
    eligible = (
        isolde.sim_manager.sim_construct.mobile_atoms if isolde.simulation_running else m.atoms
    )
    chirals = chirals[eligible.indices(chirals.chiral_atoms) != -1]
    if not len(chirals):
        raise UserError('No (mobile) chiral centres found in the selection.')
    from ..atomic.template_utils import invert_chiral_center
    # Handedness of each centre relative to its (current) target: > 0 == correct.
    oriented = numpy.sign(chirals.expected_volumes) * chirals.true_chiral_volumes

    if force:
        # Experts only: flip EVERY selected target to its mirror, then move atoms
        # only where needed. A currently-correct centre (oriented > 0) must be
        # geometry-inverted to match the flipped target; an already-inverted one
        # (oriented < 0) already matches the flipped target, so its atoms are left
        # alone -- the fewest moves to a state that minimises against all the new
        # targets.
        to_move = chirals[oriented > 0]
        sim_was_running = isolde.simulation_running
        from chimerax.isolde import session_extensions as sx
        crm = sx.get_chiral_restraint_mgr(m)
        for r in crm.add_restraints(chirals):
            r.flip_target()
        moved, unmovable = [], []
        for cc in to_move:
            (moved if invert_chiral_center(cc.chiral_atom) else unmovable).append(
                cc.chiral_atom.name
            )
        # Flip targets + geometry FIRST, then start the simulation, so its
        # automatic initial minimisation relaxes the strained result (see the
        # default-path note below).
        if not sim_was_running:
            from ..cmd import isolde_sim
            isolde_sim(session, 'start', chirals.chiral_atoms.unique_residues.atoms)
        elif moved and not isolde.sim_paused:
            # Only an actively-running (unpaused) sim overwrites the model
            # coordinates on its next step, so the inverted coords must be pushed
            # in (which also re-triggers minimisation). When paused, the coordinate
            # and (deferred) restraint-target changes are applied automatically on
            # resume, so no manual push is needed.
            isolde.sim_handler.push_coords_to_sim()
        session.logger.warning(
            'isolde chiralflip force (EXPERTS ONLY): flipped the '
            'restraint target(s) of {} chiral centre(s) to the mirror reference '
            '(inverted {} to match; {} already matched). SESSION ONLY -- these '
            'centres now read as correct in validation, but the CCD/chirals.json '
            'definition is unchanged and restored on reload, at which point these '
            'centres will be flagged as incorrect and restrained back to the '
            'canonical chirality. To persist, give the residue a new unique name -- '
            'but it will then be chirally UNRESTRAINED until you register it as a '
            'custom template.'.format(len(chirals), len(moved),
                                      len(chirals) - len(to_move))
        )
        if unmovable:
            session.logger.warning(
                'Target flipped but geometry could NOT be '
                'inverted for {} (no small rotatable substituent); these will be '
                'strained -- delete and rebuild instead.'.format(', '.join(unmovable))
            )
        return

    # Default: repair only genuinely-inverted centres by flipping their GEOMETRY
    # back to the existing target; leave correct (even if strained) ones alone.
    to_flip = chirals[oriented < 0]
    correct = chirals[oriented >= 0]
    if not len(to_flip):
        if len(chirals) == 1:
            raise UserError(
                'Chiral centre {} is already correct -- nothing to '
                'flip. (Use "force true" to deliberately invert it -- experts '
                'only.)'.format(chirals.chiral_atoms.names[0])
            )
        raise UserError(
            'All {} selected chiral centres are already correct -- '
            'nothing to flip. (Use "force true" to deliberately invert -- '
            'experts only.)'.format(len(chirals))
        )
    if len(correct):
        session.logger.warning(
            'Leaving {} already-correct chiral centre(s) '
            'untouched ({}); flipping only the {} inverted one(s).'.format(
                len(correct), ', '.join(correct.chiral_atoms.names), len(to_flip)
            )
        )
    # Flip the geometry FIRST, then start the simulation (or push into a running
    # one). A freshly-started sim runs an automatic initial energy minimisation
    # that relaxes the deliberately-strained flipped centre. Starting the sim
    # *before* flipping and pushing the inverted coords in afterwards skips that
    # relaxation -- on the very first run a post-start coordinate push does not
    # re-trigger minimisation, so dynamics resume from the strained geometry and
    # kick hard enough to corrupt neighbouring residues (the first-flip glitch).
    # Pushing into an *already-running* sim does re-trigger minimisation, so that
    # path is kept.
    sim_was_running = isolde.simulation_running
    flipped_mask = numpy.zeros(len(to_flip), dtype=bool)
    unflippable = []
    for i, cc in enumerate(to_flip):
        if invert_chiral_center(cc.chiral_atom):
            flipped_mask[i] = True
        else:
            unflippable.append(cc.chiral_atom.name)
    flipped = to_flip[flipped_mask]
    if not sim_was_running:
        from ..cmd import isolde_sim
        isolde_sim(session, 'start', to_flip.chiral_atoms.unique_residues.atoms)
    elif len(flipped) and not isolde.sim_paused:
        # Only an actively-running (unpaused) sim overwrites the model coordinates
        # on its next step, so push the inverted coords in (re-triggering
        # minimisation). When paused, the change is applied automatically on resume.
        isolde.sim_handler.push_coords_to_sim()
    if len(flipped):
        session.logger.info(
            'Flipped chiral centre(s) {}; relaxing under the '
            'chiral restraint.'.format(', '.join(flipped.chiral_atoms.names))
        )
    if unflippable:
        session.logger.warning(
            'Could not flip chiral centre(s) {} by rotation '
            '(no small rotatable substituent); delete and rebuild these instead.'.format(
                ', '.join(unflippable)
            )
        )


def register_manip_commands(logger):
    from chimerax.core.commands import (
        register,
        CmdDesc,
        BoolArg,
    )
    from chimerax.atomic import AtomsArg

    def register_isolde_pepflip():
        desc = CmdDesc(
            required=[('atoms', AtomsArg)],
            synopsis=('Attempt to flip the peptide bond N-terminal to each '
                'selected residue, starting a suitable simulation if none '
                'is currently running. Requires ISOLDE to already be '
                'initialised. Residues outside the model currently selected '
                'for ISOLDE (or outside the mobile selection if a simulation '
                'is running) will not be flipped.')
        )
        register('isolde pepflip', desc, pep_flip, logger=logger)

    def register_isolde_cisflip():
        desc = CmdDesc(
            required=[('atoms', AtomsArg)],
            synopsis=('Attempt to flip the peptide bond N-terminal to each '
                'selected residue from cis to trans or vice versa, starting a '
                'suitable simulation if none is currently running. Requires '
                'ISOLDE to already be initialised. Residues outside the model '
                'currently selected for ISOLDE (or outside the mobile selection '
                'if a simulation is running) will not be flipped.')
        )
        register('isolde cisflip', desc, cis_flip, logger=logger)

    def register_isolde_chiralflip():
        desc = CmdDesc(
            required=[('atoms', AtomsArg)],
            keyword=[('force', BoolArg)],
            synopsis=(
                'Repair inverted chiral centres by rotating the smallest '
                'substituent 180 degrees, starting a suitable simulation if none '
                'is running. Requires ISOLDE to already be initialised. By default '
                'only genuinely-inverted centres are flipped; centres outside the '
                'model (or outside the mobile selection during a simulation), and '
                'fused-ring centres with no small rotatable substituent, are not. '
                'EXPERTS ONLY: "force true" inverts every selected centre AND moves '
                'its restraint target to accept the new handedness (session-only; '
                'for building/exploring stereoisomers).'
            )
        )
        register('isolde chiralflip', desc, chiral_flip, logger=logger)


    register_isolde_pepflip()
    register_isolde_cisflip()
    register_isolde_chiralflip()
