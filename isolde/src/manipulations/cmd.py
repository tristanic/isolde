
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



def register_manip_commands(logger):
    from chimerax.core.commands import (
        register, CmdDesc,
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


    register_isolde_pepflip()
    register_isolde_cisflip()
