# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



from ..constants import defaults
OPENMM_ANGLE_UNIT = defaults.OPENMM_ANGLE_UNIT
OPENMM_RADIAL_SPRING_UNIT = defaults.OPENMM_RADIAL_SPRING_UNIT

class Peptide_Bond_Flipper:
    '''
    Responsible for flipping a peptide bond by adding temporary phi and psi
    restraints 180 degrees away from their starting angles. Automatically
    removes the restraints and deletes itself once satisfied, or when unable
    to satisfy the restraints after a number of steps defined by
    :attr:`SimParams.peptide_flipper_max_rounds`. On failure, a warning is
    printed to the ChimeraX log. If the restraints are satisfied to within
    :attr:`SimParams.dihedral_restraint_cutoff_angle` within the step limit,
    a "polishing" phase is triggered in which the restraint cutoff angles are
    reduced to zero and maintained for ten coordinate updates, at which point
    the restraints are released. There is no real need to keep a handle to a
    :class:`Peptide_Bond_Flipper` instance - just fire and forget.
    '''
    def __init__(self, isolde, residue):
        '''
        Initialisation immediately starts the peptide flipping process. Adds a
        handler to the :attr:`chimerax.AtomicStructure.triggers` 'coord update'
        trigger to check its progress and make decisions after each coordinate
        update.

        Args:
            * isolde:
                - the current :class:`Isolde` session
            * residue:
                - a currently-mobile :class:`chimerax.Residue`
        '''
        self.isolde = isolde
        if not isolde.simulation_running:
            raise TypeError('Simulation must be running!')
        session = self.session = isolde.session
        sim_params = self.sim_params = isolde.sim_params
        structure = self.structure = residue.structure
        from ..session_extensions import get_proper_dihedral_restraint_mgr
        pdr_m = self.pdr_m = get_proper_dihedral_restraint_mgr(structure)
        phi = self.phi = pdr_m.add_restraint_by_residue_and_name(residue, 'phi')
        if phi is None:
            raise TypeError('Residue does not have an N-terminal peptide bond!')
        psi = self.psi = pdr_m.add_restraint_by_residue_and_name(
            phi.dihedral.atoms[0].residue, 'psi'
        )
        from math import pi
        from ..molarray import Proper_Dihedral_Restraints
        phipsi = self.phipsi = Proper_Dihedral_Restraints([phi, psi])
        import numpy
        if numpy.any(phipsi.sim_indices == -1):
            raise TypeError('Peptide bond must be mobile in the simulation!')
        phipsi.targets = phipsi.dihedrals.angles + pi
        phipsi.spring_constants = \
            sim_params.peptide_bond_spring_constant.value_in_unit(
                OPENMM_RADIAL_SPRING_UNIT
            )
        phipsi.cutoffs = sim_params.dihedral_restraint_cutoff_angle.value_in_unit(
            OPENMM_ANGLE_UNIT
        )
        phipsi.enableds = True
        self._attempt_counter = 0
        self._coord_update_handler = isolde.sim_handler.triggers.add_handler(
            'coord update', self._initial_change_cb
        )

    def _initial_change_cb(self, *_):
        import numpy
        phipsi = self.phipsi
        isolde = self.isolde
        offsets = numpy.abs(phipsi.offsets)
        if numpy.all(offsets < phipsi.cutoffs):
            phipsi.cutoffs = 0
            self._attempt_counter = 0
            self._coord_update_handler = isolde.sim_handler.triggers.add_handler(
                'coord update', self._final_polish_cb
            )
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER
        self._attempt_counter += 1
        if self._attempt_counter > self.sim_params.peptide_flipper_max_rounds:
            self.session.logger.info(
                'Unable to flip peptide bond after {} rounds. Giving up.'\
                .format(self.sim_params.peptide_flipper_max_rounds))
            phipsi.enableds = False
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER

    def _final_polish_cb(self, *_):
        if self._attempt_counter >= 10:
            import numpy
            phipsi = self.phipsi
            phipsi.enableds = False
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER
        self._attempt_counter += 1
