# @Author: Tristan Croll
# @Date:   30-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 30-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Abstract interface for MD parameterisation backends.

:class:`ParameterisationProvider` is an ABC that any parameterisation backend
must implement.  :class:`ForceFieldParameterisationProvider` wraps the existing
XML-based OpenMM forcefield pipeline and is the only concrete implementation for
now; future backends (ML forcefields, etc.) implement the ABC without touching
:class:`SimHandler`.

Circular-import note
--------------------
This module deliberately avoids module-level imports from ``openmm_interface``.
All references to ``openmm_interface`` symbols (``find_residue_templates``,
``create_openmm_topology``, ``UnparameterisedResiduesError``) are deferred to
inside method bodies so that the ``param_provider → openmm_interface`` dependency
is a runtime-only dependency.  ``openmm_interface.py`` never imports this module.
'''

from abc import ABC, abstractmethod


class ParameterisationProvider(ABC):
    '''
    Abstract base for objects that convert a :class:`SimConstruct` +
    :class:`SimParams` into an OpenMM topology and system ready for simulation.

    The single required method is :meth:`build_system`.

    .. note::
        Implementations must **not** call ``_set_core_force_groups``.  That is
        the responsibility of :class:`SimHandler` after ``build_system`` returns,
        because the force-group constants live in ``openmm_interface.py``.
    '''

    @abstractmethod
    def build_system(self, sim_construct, sim_params, logger=None):
        '''
        Parameterise *sim_construct* and return the OpenMM objects needed to
        drive a simulation.

        Parameters
        ----------
        sim_construct : SimConstruct
            Atoms, residues, and fixed-atom bookkeeping for this run.
        sim_params : SimParams
            User-facing simulation parameters (forcefield name, cutoffs,
            constraints, etc.).
        logger : ChimeraX logger or None
            Used for warnings/info during template resolution.

        Returns
        -------
        topology : openmm.app.Topology
        system : openmm.System
            The returned system has **not** had force groups assigned;
            :class:`SimHandler` assigns them after this call returns.
        '''


class ForceFieldParameterisationProvider(ParameterisationProvider):
    '''
    Concrete provider that wraps ISOLDE's existing XML-based OpenMM forcefield
    pipeline (:class:`ForcefieldMgr` + AMBER14/CHARMM36 XML files).

    Parameters
    ----------
    forcefield_mgr : ForcefieldMgr
        Dict-like manager that lazy-loads :class:`ForceField` objects by name.
    '''

    def __init__(self, forcefield_mgr):
        self._forcefield_mgr = forcefield_mgr

    @property
    def forcefield_mgr(self):
        '''The underlying :class:`ForcefieldMgr` (for backward-compat access).'''
        return self._forcefield_mgr

    def build_system(self, sim_construct, sim_params, logger=None):
        ff = self._forcefield_mgr[sim_params.forcefield]
        ligand_db = self._forcefield_mgr.ligand_db(sim_params.forcefield)
        atoms = sim_construct.all_atoms
        all_residues = sim_construct.all_residues
        return self._build_with_template_recovery(
            ff, atoms, all_residues, sim_params, ligand_db, logger)

    def _build_with_template_recovery(self, ff, atoms, all_residues,
                                      sim_params, ligand_db, logger):
        '''
        Build the OpenMM topology and system, retrying once on stale
        ``isolde_template_name`` overrides.

        Migrated verbatim from ``SimHandler._build_system_with_template_recovery``.
        Rebuilding from scratch — rather than just re-running template matching —
        is important: ``find_residue_templates`` deliberately skips its cysteine /
        USER_ / ligand logic for residues that already carry an override, so the
        residue only gets a fair shot at automatic assignment after the override
        is gone.

        Returns ``(topology, system)``.  The system has **not** had force groups set.
        '''
        # Deferred imports to avoid a circular dependency at module load time.
        from chimerax.isolde.openmm.openmm_interface import (
            find_residue_templates,
            create_openmm_topology,
            UnparameterisedResiduesError,
        )
        for attempt in (0, 1):
            template_dict = find_residue_templates(all_residues, ff,
                ligand_db=ligand_db, logger=logger,
                clear_failed_overrides=True)
            top, residue_templates = create_openmm_topology(atoms, template_dict)
            try:
                system = self._create_system(ff, top, sim_params,
                    residue_templates, residues=all_residues)
            except UnparameterisedResiduesError as e:
                offenders = [r for r in (e.unmatched + e.ambiguous)
                    if getattr(r, 'isolde_template_name', None) is not None]
                if attempt == 0 and offenders:
                    for r in offenders:
                        if logger is not None:
                            logger.warning(
                                'Residue {} was assigned template "{}" via its '
                                'isolde_template_name override, but that template '
                                'does not match the residue. Clearing the override '
                                'and retrying with automatic template assignment.'
                                .format(r, r.isolde_template_name))
                        r.isolde_template_name = None
                    continue
                raise
            return top, system

    def _create_system(self, forcefield, top, params, residue_templates,
                       residues=None):
        '''
        Call ``forcefield.createSystem`` and return the raw :class:`openmm.System`.

        Migrated from ``SimHandler._create_openmm_system`` with one intentional
        change: the ``_set_core_force_groups`` call is removed.  Force-group
        assignment is :class:`SimHandler`'s responsibility.
        '''
        from chimerax.isolde.openmm.openmm_interface import UnparameterisedResiduesError

        residue_to_template, ambiguous, unassigned = forcefield.assignTemplates(
            top, ignoreExternalBonds=True, explicit_templates=residue_templates
        )
        if len(ambiguous) or len(unassigned):
            # Map offending OpenMM topology residues back to ChimeraX residues
            # (OpenMM res.index matches the order added in create_openmm_topology,
            # which matches `residues`).
            unmatched_cx = []
            ambiguous_cx = []
            if residues is not None:
                unmatched_cx = [residues[r.index] for r in unassigned]
                ambiguous_cx = [residues[r.index] for r in ambiguous]
            raise UnparameterisedResiduesError(unmatched=unmatched_cx,
                ambiguous=ambiguous_cx)

        system_params = {
            'nonbondedMethod':      params.nonbonded_cutoff_method,
            'nonbondedCutoff':      params.nonbonded_cutoff_distance,
            'constraints':          params.rigid_bonds,
            'rigidWater':           params.rigid_water,
            'removeCMMotion':       params.remove_c_of_m_motion,
            'ignoreExternalBonds':  True,
            'residueTemplates':     residue_templates,
        }
        return forcefield.createSystem(top, **system_params)
