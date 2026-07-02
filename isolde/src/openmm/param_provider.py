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


def _garnet_available() -> bool:
    '''True when garnetff is importable.'''
    try:
        import garnetff  # noqa: F401
        return True
    except ImportError:
        return False


def _espaloma_available() -> bool:
    '''True when espaloma, rdkit, and openff-toolkit are all importable.'''
    try:
        import espaloma  # noqa: F401
        from rdkit import Chem  # noqa: F401
        from openff.toolkit import Molecule  # noqa: F401
        return True
    except ImportError:
        return False


def _openff_available() -> bool:
    '''True when openff-toolkit and rdkit are importable.'''
    try:
        from openff.toolkit import ForceField, Molecule  # noqa: F401
        from rdkit import Chem  # noqa: F401
        return True
    except ImportError:
        return False


def _mmff_available() -> bool:
    '''True when rdkit is importable (MMFF94 is bundled with rdkit).'''
    try:
        from rdkit.Chem import rdForceFieldHelpers  # noqa: F401
        return True
    except ImportError:
        return False


def available_parameterisation_backends() -> list:
    '''
    Return the list of backend names that can be passed to
    ``sim_params.forcefield``.
    '''
    backends = ['amber14']
    if _garnet_available():
        backends.append('amber14+garnet')
    if _espaloma_available():
        backends.append('amber14+espaloma')
    if _openff_available():
        backends.append('amber14+openff')
    if _mmff_available():
        backends.append('amber14+mmff')
    backends.append('charmm36')
    return backends


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


def _residue_structural_signature(residue):
    '''
    A string fingerprint of a ChimeraX residue's atom names/elements and
    internal bonding, order-independent. Used (alongside a forcefield
    identity key and any explicit template-name override) to detect whether
    a residue's cached AMBER type assignment is still valid -- catches
    silent atom/bond edits that no ``find_residue_templates`` heuristic
    reacts to. Deliberately not hashed: Python's ``hash()`` of strings is
    process-salted, so a hash computed in one process would not reliably
    compare equal after a restart or session reload.
    '''
    atoms = residue.atoms
    atom_part = tuple(sorted(zip(atoms.names.tolist(), atoms.element_names.tolist())))
    b0, b1 = atoms.intra_bonds.atoms
    bond_part = tuple(sorted(
        tuple(sorted((n1, n2))) for n1, n2 in zip(b0.names.tolist(), b1.names.tolist())
    ))
    return repr((atom_part, bond_part))


import weakref
_amber_cache_attrs_registered_sessions = weakref.WeakSet()


def _ensure_amber_cache_attrs_registered(session):
    '''
    ``Atom.register_attr``/``Residue.register_attr`` calls normally happen via
    the bundle's ``custom-init`` hook (see ``chimerax.isolde.register_amber_type_cache_attrs``),
    but that hook is only guaranteed to fire during a normal interactive
    ChimeraX startup -- it does not reliably run before a headless
    ``--script``/``runscript`` invocation (confirmed empirically: reading an
    unregistered custom attribute raises ``AttributeError`` even after opening
    a structure in such a session). Guard against that here rather than
    silently miscaching, since ISOLDE is also driven headlessly via its
    MCP/agent control surface.

    Tracked via a ``WeakSet`` rather than ``id(session)`` in a plain ``set``:
    a plain set of ids would never release entries for destroyed sessions, and
    -- since ``id()`` is just a memory address -- a later, unrelated session
    object could be allocated at the same address and be wrongly treated as
    already registered.
    '''
    if session in _amber_cache_attrs_registered_sessions:
        return
    from chimerax.isolde import register_amber_type_cache_attrs
    register_amber_type_cache_attrs(session)
    _amber_cache_attrs_registered_sessions.add(session)


def _forcefield_identity_key(forcefield, ff_name):
    '''
    A string identifying "this forcefield, as currently loaded" -- changes
    whenever the forcefield is switched, ISOLDE/OpenMM is upgraded (mirroring
    :class:`ForcefieldMgr`'s own pickle-cache invalidation check), or a new
    residue template is registered at runtime (glycan/ligand-db/GARNET/user
    ``loadFile`` calls all bump :attr:`ForceField._isolde_template_generation`
    via the overridden ``registerResidueTemplate``).
    '''
    from openmm import version as openmm_version
    from chimerax.isolde import __version__ as isolde_version
    generation = getattr(forcefield, '_isolde_template_generation', 0)
    return '{}|{}|{}|{}'.format(ff_name, openmm_version.version, isolde_version, generation)


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
                    residue_templates, residues=all_residues, atoms=atoms)
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
                       residues, atoms):
        '''
        Build the OpenMM :class:`System` for *top*, using cached per-atom AMBER
        types where available and only running full OpenMM template matching
        (the expensive, isomorphism-based step) for residues whose cache is
        missing or stale.

        ``residues``/``atoms`` are the ChimeraX ``Residues``/``Atoms`` arrays
        used to build *top* via ``create_openmm_topology``. Their ordering is
        guaranteed to align positionally with ``top.residues()``/``top.atoms()``
        -- ``residues[i]``/``atoms[i]`` is the ChimeraX counterpart of OpenMM
        residue/atom index ``i`` (see ``SimConstruct`` and
        ``create_openmm_topology``; this is the same invariant the old
        unmatched/ambiguous-residue error mapping below already relied on).

        Migrated from ``SimHandler._create_openmm_system``. Force-group
        assignment remains :class:`SimHandler`'s responsibility.
        '''
        from chimerax.isolde.openmm.openmm_interface import UnparameterisedResiduesError
        import json

        if len(residues):
            _ensure_amber_cache_attrs_registered(residues.unique_structures[0].session)

        ff_key = _forcefield_identity_key(forcefield, params.forcefield)
        top_residues = list(top.residues())
        data = forcefield._SystemData(top)

        miss = []  # list of (omm_res, cx_res, explicit_name)
        n_hit = 0
        for cx_res, omm_res in zip(residues, top_residues):
            explicit_name = residue_templates.get(omm_res, None)
            expected_key = '{}|{}|{}'.format(
                ff_key, explicit_name, _residue_structural_signature(cx_res))
            omm_res_atoms = list(omm_res.atoms())
            cx_res_atoms = [atoms[a.index] for a in omm_res_atoms]
            # Direct attribute access raises AttributeError for a *registered*
            # custom attribute that has simply never been set on this instance
            # -- registration declares the attribute, it does not give unset
            # instances a default value. getattr(..., None) is the existing
            # codebase convention for exactly this reason (cf. find_residue_templates'
            # use of isolde_template_name).
            if (getattr(cx_res, 'isolde_amber_cache_key', None) == expected_key
                    and all(getattr(a, 'isolde_amber_type', None) is not None for a in cx_res_atoms)):
                for omm_atom, cx_atom in zip(omm_res_atoms, cx_res_atoms):
                    data.atomType[omm_atom] = cx_atom.isolde_amber_type
                    p = getattr(cx_atom, 'isolde_amber_params', None)
                    data.atomParameters[omm_atom] = json.loads(p) if p else {}
                n_hit += 1
            else:
                miss.append((omm_res, cx_res, explicit_name))

        # Recorded before attempting to resolve the misses, so that a
        # residue which turns out to be genuinely unparameterisable (e.g. an
        # atom deletion left it with no matching template) still reports
        # accurate hit/miss counts to the caller/tests, rather than stale
        # numbers from a previous successful build.
        self.last_build_stats = {
            'residues_total': len(top_residues),
            'residues_cache_hit': n_hit,
            'residues_cache_miss': len(miss),
        }

        if miss:
            miss_top_residues = [m[0] for m in miss]
            residue_to_template, ambiguous, unassigned = forcefield.assignTemplates(
                top, ignoreExternalBonds=True, explicit_templates=residue_templates,
                residues=miss_top_residues
            )
            if len(ambiguous) or len(unassigned):
                # Map offending OpenMM topology residues back to ChimeraX residues
                # (OpenMM res.index matches the order added in create_openmm_topology,
                # which matches `residues`).
                unmatched_cx = [residues[r.index] for r in unassigned]
                ambiguous_cx = [residues[r.index] for r in ambiguous]
                raise UnparameterisedResiduesError(unmatched=unmatched_cx,
                    ambiguous=ambiguous_cx)

            cx_res_and_name_by_top_res = {m[0]: (m[1], m[2]) for m in miss}
            for omm_res, (template, match) in residue_to_template.items():
                data.recordMatchedAtomParameters(omm_res, template, match)
                cx_res, explicit_name = cx_res_and_name_by_top_res[omm_res]
                if template.virtualSites:
                    # None of ISOLDE's currently-loaded AMBER14 XML files define
                    # virtual sites, but if one ever does, never cache it --
                    # always re-resolve rather than risk a stale/incomplete
                    # cache entry silently producing a wrong system. The
                    # current build is still correct either way; only future
                    # cache reuse is affected.
                    cx_res.isolde_amber_cache_key = None
                    continue
                cache_key = '{}|{}|{}'.format(
                    ff_key, explicit_name, _residue_structural_signature(cx_res))
                omm_res_atoms = list(omm_res.atoms())
                cx_res_atoms = [atoms[a.index] for a in omm_res_atoms]
                for omm_atom, cx_atom in zip(omm_res_atoms, cx_res_atoms):
                    cx_atom.isolde_amber_type = data.atomType[omm_atom]
                    p = data.atomParameters.get(omm_atom) or {}
                    cx_atom.isolde_amber_params = json.dumps(p) if p else None
                cx_res.isolde_amber_cache_key = cache_key

        return self._build_system_from_data(forcefield, data, top, params)

    def _build_system_from_data(self, forcefield, data, top, params):
        '''
        Construct the OpenMM :class:`System` from a fully-populated
        ``_SystemData`` (``atomType``/``atomParameters`` set for every atom,
        either from cache or from fresh template matching in
        :meth:`_create_system`).

        This is a close adaptation of the non-matching "tail" of
        ``openmm.app.forcefield.ForceField.createSystem`` (vendored under
        ``extern/openmm/wrappers/python/openmm/app/forcefield.py``, roughly
        lines 1267-1423 as of the OpenMM version pinned by ``ForcefieldMgr``).
        It is copied close to verbatim -- rather than re-derived -- because the
        angle/proper/improper enumeration order affects floating-point
        summation order in the resulting energy, and because the actual AMBER
        math lives entirely in the force-generator objects
        (``forcefield._forces``), which are called completely unmodified here.
        If ``ForcefieldMgr._openmm_version`` is ever bumped to a version with a
        changed ``createSystem``, re-diff this method against the new tail.
        '''
        import openmm as mm
        from openmm import unit
        from openmm.app import element as elem
        from openmm.app.forcefield import NoCutoff, CutoffNonPeriodic, AllBonds, HAngles, HBonds
        import itertools

        nonbondedMethod = params.nonbonded_cutoff_method
        nonbondedCutoff = params.nonbonded_cutoff_distance
        constraints = params.rigid_bonds
        rigidWater = params.rigid_water
        removeCMMotion = params.remove_c_of_m_motion
        args = {'switchDistance': None, 'flexibleConstraints': False,
                'drudeMass': 0.4 * unit.amu}

        if rigidWater is None:
            # ISOLDE's SimParams always supplies a concrete bool; this would
            # only trigger via direct/test misuse of this method.
            raise ValueError('_build_system_from_data requires an explicit '
                'rigid_water setting (True/False), not None.')
        rigidResidue = [False] * top.getNumResidues()
        if rigidWater:
            for res in top.residues():
                if res.name == 'HOH':
                    rigidResidue[res.index] = True

        # Calculate atom classes once for reuse by all generators (important
        # for Amoeba-style forces). Present on the currently-installed OpenMM
        # but not on the vendored extern/openmm reference copy this method is
        # otherwise adapted from -- a reminder that the installed runtime
        # version, not the vendored copy, is the actual source of truth to
        # re-diff against.
        data.setAtomClasses(forcefield)

        sys = mm.System()
        for atom in top.atoms():
            if atom not in data.atomType:
                raise Exception("Could not identify atom type for atom '%s'." % str(atom))
            typename = data.atomType[atom]
            if typename not in forcefield._atomTypes:
                msg = "Could not find typename '%s' for atom '%s' in list of known atom types.\n" % (typename, str(atom))
                msg += "Known atom types are: %s" % str(forcefield._atomTypes.keys())
                raise Exception(msg)
            mass = forcefield._atomTypes[typename].mass
            sys.addParticle(mass)

        boxVectors = top.getPeriodicBoxVectors()
        if boxVectors is not None:
            sys.setDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2])
        elif nonbondedMethod not in (NoCutoff, CutoffNonPeriodic):
            raise ValueError('Requested periodic boundary conditions for a Topology '
                              'that does not specify periodic box dimensions')

        # Make a list of all unique angles

        uniqueAngles = set()
        for bond in data.bonds:
            for atom in data.bondedToAtom[bond.atom1]:
                if atom != bond.atom2:
                    if atom < bond.atom2:
                        uniqueAngles.add((atom, bond.atom1, bond.atom2))
                    else:
                        uniqueAngles.add((bond.atom2, bond.atom1, atom))
            for atom in data.bondedToAtom[bond.atom2]:
                if atom != bond.atom1:
                    if atom > bond.atom1:
                        uniqueAngles.add((bond.atom1, bond.atom2, atom))
                    else:
                        uniqueAngles.add((atom, bond.atom2, bond.atom1))
        data.angles = sorted(list(uniqueAngles))

        # Make a list of all unique proper torsions

        uniquePropers = set()
        for angle in data.angles:
            for atom in data.bondedToAtom[angle[0]]:
                if atom not in angle:
                    if atom < angle[2]:
                        uniquePropers.add((atom, angle[0], angle[1], angle[2]))
                    else:
                        uniquePropers.add((angle[2], angle[1], angle[0], atom))
            for atom in data.bondedToAtom[angle[2]]:
                if atom not in angle:
                    if atom > angle[0]:
                        uniquePropers.add((angle[0], angle[1], angle[2], atom))
                    else:
                        uniquePropers.add((atom, angle[2], angle[1], angle[0]))
        data.propers = sorted(list(uniquePropers))

        # Make a list of all unique improper torsions

        for atom in range(len(data.bondedToAtom)):
            bondedTo = data.bondedToAtom[atom]
            if len(bondedTo) > 2:
                for subset in itertools.combinations(bondedTo, 3):
                    data.impropers.append((atom, subset[0], subset[1], subset[2]))

        # Identify bonds that should be implemented with constraints

        if constraints == AllBonds or constraints == HAngles:
            for bond in data.bonds:
                bond.isConstrained = True
        elif constraints == HBonds:
            for bond in data.bonds:
                atom1 = data.atoms[bond.atom1]
                atom2 = data.atoms[bond.atom2]
                bond.isConstrained = atom1.element is elem.hydrogen or atom2.element is elem.hydrogen
        for bond in data.bonds:
            atom1 = data.atoms[bond.atom1]
            atom2 = data.atoms[bond.atom2]
            if rigidResidue[atom1.residue.index] and rigidResidue[atom2.residue.index]:
                bond.isConstrained = True

        # Identify angles that should be implemented with constraints

        if constraints == HAngles:
            for angle in data.angles:
                atom1 = data.atoms[angle[0]]
                atom2 = data.atoms[angle[1]]
                atom3 = data.atoms[angle[2]]
                numH = 0
                if atom1.element is elem.hydrogen:
                    numH += 1
                if atom3.element is elem.hydrogen:
                    numH += 1
                data.isAngleConstrained.append(numH == 2 or (numH == 1 and atom2.element is elem.oxygen))
        else:
            data.isAngleConstrained = len(data.angles) * [False]
        for i in range(len(data.angles)):
            angle = data.angles[i]
            atom1 = data.atoms[angle[0]]
            atom2 = data.atoms[angle[1]]
            atom3 = data.atoms[angle[2]]
            if (rigidResidue[atom1.residue.index] and rigidResidue[atom2.residue.index]
                    and rigidResidue[atom3.residue.index]):
                data.isAngleConstrained[i] = True

        # Add virtual sites

        for atom in data.virtualSites:
            (site, site_atoms, excludeWith) = data.virtualSites[atom]
            index = atom.index
            data.excludeAtomWith[excludeWith].append(index)
            if site.type == 'average2':
                sys.setVirtualSite(index, mm.TwoParticleAverageSite(site_atoms[0], site_atoms[1], site.weights[0], site.weights[1]))
            elif site.type == 'average3':
                sys.setVirtualSite(index, mm.ThreeParticleAverageSite(site_atoms[0], site_atoms[1], site_atoms[2], site.weights[0], site.weights[1], site.weights[2]))
            elif site.type == 'outOfPlane':
                sys.setVirtualSite(index, mm.OutOfPlaneSite(site_atoms[0], site_atoms[1], site_atoms[2], site.weights[0], site.weights[1], site.weights[2]))
            elif site.type == 'localCoords':
                sys.setVirtualSite(index, mm.LocalCoordinatesSite(site_atoms, site.originWeights, site.xWeights, site.yWeights, site.localPos))

        # Add forces to the System

        for force in forcefield._forces:
            force.createForce(sys, data, nonbondedMethod, nonbondedCutoff, args)
        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())

        # Let force generators do postprocessing

        for force in forcefield._forces:
            if hasattr(force, 'postprocessSystem'):
                force.postprocessSystem(sys, data, args)

        # Execute scripts found in the XML files.

        for script in forcefield._scripts:
            exec(script, locals())

        return sys


# ---------------------------------------------------------------------------
# Base class for AMBER14 + external-backend providers
# ---------------------------------------------------------------------------

class LigandBackedParameterisationProvider(ForceFieldParameterisationProvider):
    '''
    Base for providers that combine AMBER14 for standard residues with an
    external backend (GARNET, Espaloma, OpenFF, …) for non-polymer ligands.

    Subclasses must implement three hooks:

    - :attr:`backend_forcefield_key` — the ``sim_params.forcefield`` string
      this provider owns (e.g. ``'amber14+garnet'``).
    - :meth:`_residue_is_eligible` — whether to send a residue to the
      backend.  Default: any non-polymer residue with ≥ 3 heavy atoms.
    - :meth:`_parameterise_one` — run the backend on one residue, returning
      ``(ffxml_str, extra_dict)``.  *extra_dict* is merged into the cache
      entry and made available to :meth:`_patch_system`.

    Subclasses may override :meth:`_patch_system` (default: no-op) to
    post-process the assembled system (e.g. replace LJ with a custom VdW).
    '''

    def __init__(self, forcefield_mgr):
        super().__init__(forcefield_mgr)
        # {resname: {'xml': str, **backend_extra}}
        self._template_cache = {}

    # -- abstract interface --

    @property
    def backend_forcefield_key(self):
        raise NotImplementedError

    def _residue_is_eligible(self, residue):
        from chimerax.atomic import Residue
        if residue.polymer_type != Residue.PT_NONE:
            return False
        heavy = [a for a in residue.atoms if a.element.number != 1]
        return len(heavy) >= 3

    def _parameterise_one(self, residue, logger):
        raise NotImplementedError

    def _patch_system(self, system, top):
        pass

    # -- routing --

    def build_system(self, sim_construct, sim_params, logger=None):
        if sim_params.forcefield != self.backend_forcefield_key:
            return super().build_system(sim_construct, sim_params, logger)
        ff  = self._forcefield_mgr['amber14']
        ldb = self._forcefield_mgr.ligand_db('amber14')
        return self._build_with_ligand_fallback(
            ff, sim_construct.all_atoms, sim_construct.all_residues,
            sim_params, ldb, logger)

    # -- shared machinery --

    def _build_with_ligand_fallback(self, ff, atoms, all_residues,
                                    sim_params, ligand_db, logger):
        '''
        3-attempt retry loop.

        0.  Stale ``isolde_template_name`` recovery.
        1.  Backend fill-in for residues still unparameterised after
            :meth:`_pre_parameterise_ligands`.
        2.  Propagate error.
        '''
        self._pre_parameterise_ligands(all_residues, ff, logger)
        from chimerax.isolde.openmm.openmm_interface import (
            find_residue_templates,
            create_openmm_topology,
            UnparameterisedResiduesError,
        )
        for attempt in (0, 1, 2):
            template_dict = find_residue_templates(
                all_residues, ff,
                ligand_db=ligand_db,
                logger=logger,
                clear_failed_overrides=True)
            top, residue_templates = create_openmm_topology(atoms, template_dict)
            try:
                system = self._create_system(
                    ff, top, sim_params, residue_templates,
                    residues=all_residues, atoms=atoms)
            except UnparameterisedResiduesError as e:
                stale = [r for r in (e.unmatched + e.ambiguous)
                         if getattr(r, 'isolde_template_name', None) is not None]
                if attempt == 0 and stale:
                    for r in stale:
                        if logger is not None:
                            logger.warning(
                                'Residue {} was assigned template "{}" via its '
                                'isolde_template_name override, but that template '
                                'does not match the residue. Clearing the override '
                                'and retrying with automatic template assignment.'
                                .format(r, r.isolde_template_name))
                        r.isolde_template_name = None
                    continue
                if attempt == 1 and e.unmatched:
                    self._parameterise_unmatched(e.unmatched, ff, logger)
                    continue
                raise
            self._patch_system(system, top)
            return top, system

    def _pre_parameterise_ligands(self, all_residues, forcefield, logger):
        '''
        Proactively parameterise all eligible non-polymer residues and register
        their templates under the ``USER_`` prefix so they take precedence over
        any AMBER14 core-XML template during the subsequent
        ``find_residue_templates`` call.
        '''
        import io
        backend = self.backend_forcefield_key.split('+')[1].upper()

        for r in all_residues:
            if not self._residue_is_eligible(r):
                continue

            name = r.name
            existing = getattr(r, 'isolde_template_name', None)
            if existing is not None and not existing.startswith('USER_'):
                r.isolde_template_name = None

            if name not in self._template_cache:
                try:
                    if logger is not None:
                        logger.status(f'{backend}: parameterising {name}...')
                    xml, extra = self._parameterise_one(r, logger)
                    self._template_cache[name] = {'xml': xml, **extra}
                    if logger is not None:
                        logger.info(f'{backend}: parameterisation complete for {name}.')
                except Exception as exc:
                    if logger is not None:
                        logger.warning(
                            f'{backend}: could not parameterise {name} ({exc}); '
                            f'falling back to AMBER14.')
                    continue

            forcefield.loadFile(
                io.StringIO(self._template_cache[name]['xml']),
                resname_prefix='USER_')

    def _parameterise_unmatched(self, unmatched_residues, forcefield, logger):
        '''
        Parameterise still-unmatched non-polymer residues on the second retry.
        '''
        import io
        from chimerax.atomic import Residue
        backend = self.backend_forcefield_key.split('+')[1].upper()

        for r in unmatched_residues:
            if r.polymer_type != Residue.PT_NONE:
                continue
            name = r.name
            if name not in self._template_cache:
                if logger is not None:
                    logger.status(f'{backend}: parameterising {name}...')
                xml, extra = self._parameterise_one(r, logger)
                self._template_cache[name] = {'xml': xml, **extra}
                if logger is not None:
                    logger.info(f'{backend}: parameterisation complete for {name}.')
            forcefield.loadFile(
                io.StringIO(self._template_cache[name]['xml']),
                resname_prefix='USER_')


# ---------------------------------------------------------------------------
# Dispatcher — routes to whichever backend owns the current forcefield key
# ---------------------------------------------------------------------------

class _DispatchingProvider(ForceFieldParameterisationProvider):
    '''
    Routes each ``sim_params.forcefield`` value to the appropriate
    :class:`LigandBackedParameterisationProvider`.  Created by
    :func:`make_parameterisation_provider` when multiple backends are installed.
    '''

    def __init__(self, forcefield_mgr, backend_providers):
        super().__init__(forcefield_mgr)
        self._backends = {p.backend_forcefield_key: p
                          for p in backend_providers}

    def build_system(self, sim_construct, sim_params, logger=None):
        backend = self._backends.get(sim_params.forcefield)
        if backend is not None:
            return backend.build_system(sim_construct, sim_params, logger)
        return super().build_system(sim_construct, sim_params, logger)


def make_parameterisation_provider(forcefield_mgr):
    '''
    Factory — returns the most capable parameterisation provider available.

    * No ML backends installed → :class:`ForceFieldParameterisationProvider`.
    * One ML backend installed → that backend's provider directly.
    * Multiple ML backends installed → :class:`_DispatchingProvider`.
    '''
    backends = []
    if _garnet_available():
        from .garnet_provider import GarnetParameterisationProvider
        backends.append(GarnetParameterisationProvider(forcefield_mgr))
    if _espaloma_available():
        from .espaloma_provider import EspalomaParameterisationProvider
        backends.append(EspalomaParameterisationProvider(forcefield_mgr))
    if _openff_available():
        from .openff_provider import OpenFFParameterisationProvider
        backends.append(OpenFFParameterisationProvider(forcefield_mgr))
    if _mmff_available():
        from .mmff_provider import MMFFParameterisationProvider
        backends.append(MMFFParameterisationProvider(forcefield_mgr))

    if not backends:
        return ForceFieldParameterisationProvider(forcefield_mgr)
    if len(backends) == 1:
        return backends[0]
    return _DispatchingProvider(forcefield_mgr, backends)
