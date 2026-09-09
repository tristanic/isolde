# @Author: Tristan Croll
# @Date:   24-Aug-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 24-Aug-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Whole-model garnet force-field parameters, cached keyed to ChimeraX atoms.

garnet parameterises from topology alone, so we run it **once** over the
complete model and store every term keyed to the ChimeraX :class:`Atom`
objects (per-atom charge/sigma/epsilon; per-bond/-angle/-torsion/-improper
parameters keyed by atom tuples). :meth:`GarnetParameters.subset_for` then
projects that onto any simulation's ``all_atoms`` array, remapping to OpenMM
particle indices and dropping terms that reach outside the selection.

Keying by *atom identity* (not featuriser index) is what makes
"parameterise the whole model, simulate a subset" correct: the fixed buffer's
chemically-incomplete outer edge is parameterised in full chemical context, and
carving out the simulation is a pure lookup.

This is the interim Python analogue of the C++ ``FFAtom``/``FFBond`` store
described in ``FORCEFIELD_ON_MODEL_BRIEFING.md``.
'''
import os
import numpy


# Newest committed garnet weights (WIP DTR-SF; 38-wide features). Same architecture
# as the earlier r5c round -- a drop-in retrain (r5d adds crystal-pressure/thermal
# training terms; the inference featurize/predict path is unchanged). Not wired to
# any default-pointer constant in garnet_core, so we name it and let garnet resolve
# it. Override with $ISOLDE_GARNET_CHECKPOINT or by passing checkpoint_path
# explicitly (later rounds/epochs land as dtr_sf_<round>_ep{N}.pt).
#
# A BARE FILENAME, not a path: where the weights live is garnet's business, and it
# differs between an installed wheel (garnet_core/trained_models/) and a repo
# checkout (a sibling garnetff/trained_models/). ISOLDE used to hard-code the
# checkout layout and walk up out of the garnet_core package to find it, which held
# only while the two directories happened to be siblings on disk.
_DEFAULT_CHECKPOINT_NAME = 'dtr_sf_r5d_ep1.pt'


_GARNET_MISSING_MSG = (
    "The GARNET force field requires the 'garnet-isolde' package, which is not installed "
    "in ChimeraX's Python.\n"
    "\n"
    "Install the wheel from ChimeraX's command line:\n"
    "    pip install /path/to/garnet_isolde-<version>-py3-none-any.whl\n"
    "\n"
    "Its torch / torch_geometric dependencies are pulled in automatically. To carry on "
    "without GARNET, choose one of the AMBER force fields in ISOLDE's General tab."
)


def require_garnet_core():
    '''
    Import and return ``garnet_core``, or raise with instructions.

    The bundle deliberately imports garnet lazily, so the failure mode for a missing
    install is an ImportError raised at simulation start -- far from the force-field
    selection that caused it. This turns that into something actionable.
    '''
    try:
        import garnet_core
    except ImportError as e:
        raise ImportError(_GARNET_MISSING_MSG) from e
    return garnet_core


def _resolve_bundled(name):
    '''
    Absolute path of a checkpoint shipped with garnet, looked up by bare filename.

    Delegates to ``garnet_core.weights``, which knows both the installed and the
    checkout layout. Falls back to the historical sibling walk for a ``garnet_core``
    predating that module, so an older checkout on another machine still works.
    '''
    garnet_core = require_garnet_core()
    try:
        from garnet_core import weights
    except ImportError:
        root = os.path.dirname(os.path.dirname(os.path.abspath(garnet_core.__file__)))
        return os.path.join(root, 'garnetff', 'trained_models', os.path.basename(name))
    return weights.resolve(name)


def default_checkpoint_path():
    '''Absolute path to the default garnet checkpoint (env override honoured).'''
    override = os.environ.get('ISOLDE_GARNET_CHECKPOINT')
    if override:
        return override
    return _resolve_bundled(_DEFAULT_CHECKPOINT_NAME)


def resolve_checkpoint_path(checkpoint_path=None):
    '''
    Resolve a force-field handle's ``checkpoint_path`` to an absolute path.

    * ``None`` -> the default checkpoint (``$ISOLDE_GARNET_CHECKPOINT`` honoured);
      this is the bare ``garnet`` alias's behaviour, unchanged from before.
    * an absolute path, or a relative path that already exists from the cwd ->
      used verbatim (explicit user/dev override).
    * a bare checkpoint filename (how the force-field registry pins each
      ``garnet-{run}`` variant, e.g. ``dtr_sf_r10b_ep2.pt``) -> resolved among the
      checkpoints garnet ships, the same way ``default_checkpoint_path`` resolves
      the default. A legacy repo-relative pin
      (``garnetff/trained_models/dtr_sf_r10b_ep2.pt``) also resolves, since only the
      basename is used. An explicit variant checkpoint deliberately does NOT consult
      the env override, so selecting ``garnet-r10b`` is deterministic.
    '''
    if checkpoint_path is None:
        return default_checkpoint_path()
    if os.path.isabs(checkpoint_path) or os.path.exists(checkpoint_path):
        return checkpoint_path
    return _resolve_bundled(checkpoint_path)


def _to_numpy(x):
    '''torch tensor / array-like -> detached float numpy array.'''
    if hasattr(x, 'detach'):
        x = x.detach().cpu().numpy()
    return numpy.asarray(x, dtype=float)


class SubsetParameters:
    '''
    garnet parameters projected onto one simulation's ``all_atoms`` order.

    All index fields are OpenMM particle indices (= position in ``all_atoms``).
    ``bonds``/``angles``/``propers``/``impropers`` are lists of ``(indices...,
    params...)`` tuples; per-atom arrays are numpy, length ``n_particles``.
    '''
    def __init__(self, n_particles, charges, sigmas, epsilons,
                 bonds, angles, propers, impropers, has_oop, globals_, missing_atoms,
                 bee=None, cg_b=None):
        self.n_particles = n_particles
        self.charges = charges
        self.sigmas = sigmas
        self.epsilons = epsilons
        self.bonds = bonds            # [(i, j, k, r0), ...]
        self.angles = angles          # [(i, j, k, k_theta, theta0), ...]
        self.propers = propers        # [(i, j, k, l, [k1..k6]), ...]
        self.impropers = impropers    # [(kind, (i, j, k, l), payload), ...]
        self.has_oop = has_oop        # True -> impropers are harmonic OOP, else periodic
        self.globals = globals_       # {'coulomb14scale', 'dexp_beta', ...[+form-specific]}
        self.missing_atoms = missing_atoms   # atoms with no cached parameters (should be empty)
        # Optional per-particle arrays (length n_particles) for later garnet forms, or None
        # when the feature is off for this checkpoint. ``bee``: per-atom repulsive-wall decay
        # (per-atom-wall form). ``cg_b``: per-atom Coulomb-guard decay (guard form).
        self.bee = bee
        self.cg_b = cg_b


class GarnetParameters:
    '''
    Cache of garnet FF parameters for one ChimeraX model, keyed to its atoms.

    Call :meth:`parameterise` once (or use :func:`get_garnet_parameters`, which
    caches on the model), then :meth:`subset_for` per simulation.
    '''
    def __init__(self, structure, checkpoint_path=None):
        self.structure = structure
        self.checkpoint_path = resolve_checkpoint_path(checkpoint_path)
        self._per_atom = {}       # Atom -> (charge, sigma, epsilon)
        # Optional per-atom quantities that only later garnet incarnations predict
        # (empty dict == the feature is off for this checkpoint; the builder keys off
        # emptiness). ``_dexp_b``: per-atom repulsive-wall decay for the per-atom-wall
        # form (``use_b_head``). ``_cg_b``: per-atom short-range Coulomb-guard decay
        # (``use_coulomb_guard``). See garnet_core.openmm_build._per_atom_wall / _coulomb_guard.
        self._dexp_b = {}         # Atom -> float
        self._cg_b = {}           # Atom -> float
        self._bonds = {}          # frozenset({a, b}) -> (k, r0)
        self._angles = {}         # (a_i, a_j, a_k) -> (k_theta, theta0)   (a_j = vertex)
        self._propers = {}        # (a_i, a_j, a_k, a_l) -> [k1..k6]
        self._impropers = {}      # (a_i, a_j, a_k, a_l) -> ('oop', k) | ('periodic', [k1, k2])
        self._globals = {}
        self.has_oop = False
        self.warnings = []
        self.parameterised = False

    # ------------------------------------------------------------------ build
    def parameterise(self, logger=None):
        '''Run garnet over the whole model and populate the cache.'''
        from .primitives import extract_primitives
        primitives, atoms = extract_primitives(self.structure)

        # Lazy imports: garnet_core / torch need only be present at sim time.
        # require_garnet_core() first, so a missing install reports what to install
        # instead of an unqualified "No module named 'garnet_core'".
        require_garnet_core()
        import torch
        from garnet_core.featurize import featurize
        from garnet_core.idatm_charges import infer_formal_charges
        from garnet_core.params import predict_parameters
        from garnet_core.checkpoint import load_any

        Z = primitives['atomic_numbers']
        idt = primitives['idatm_types']
        bonds = [tuple(b) for b in primitives['bonds']]
        rings = primitives.get('aromatic_rings')

        fc, warns = infer_formal_charges(Z, idt, bonds, aromatic_rings=rings)
        est = primitives.get('estimate_net_charge')
        net = int(sum(fc))
        if isinstance(est, (int, float)) and net != est:
            warns.append('net formal charge {} != estimate_net_charge {} (delta {}); '
                         'check for incomplete or mistyped residues'.format(net, est, net - est))
        self.warnings = warns
        if logger is not None:
            logger.info('garnet: parameterising {} atoms, net formal charge {}'.format(len(Z), net))
            for w in warns:
                logger.warning('garnet: {}'.format(w))

        # 38-wide features (ring + planar) are REQUIRED by the DTR-SF checkpoint.
        mol = {'atomic_numbers': Z, 'formal_charges': fc, 'aromatic': [0] * len(Z),
               'bonds': bonds, 'idatm_types': idt}
        data = featurize([mol], ring_features=True, planar_feature=True)

        model = load_any(self.checkpoint_path)
        model.eval()
        with torch.no_grad():
            params = predict_parameters(model, data)

        self._store(data, params, atoms)
        self.parameterised = True
        return self

    def _store(self, data, params, atoms):
        # For a single molecule dict, featurise index i == atoms[i] (order preserved),
        # so featuriser indices map straight back to ChimeraX atoms.
        q = _to_numpy(params['partial_charges'])
        sig = _to_numpy(params['sigma'])
        eps = _to_numpy(params['epsilon'])
        for i, a in enumerate(atoms):
            self._per_atom[a] = (float(q[i]), float(sig[i]), float(eps[i]))

        # Per-atom repulsive-wall decay (per-atom-wall form only). Present iff the
        # checkpoint's `use_b_head` is live; earlier rounds use the scalar `dexp_alpha`
        # global instead (stored below) and omit this key entirely.
        if 'dexp_b' in params:
            bee = _to_numpy(params['dexp_b'])
            for i, a in enumerate(atoms):
                self._dexp_b[a] = float(bee[i])
        # Per-atom short-range Coulomb-guard decay (guard form only). garnet_core emits
        # `coulomb_guard_b` (per-atom) + `coulomb_guard_w` (scalar) together, and only when
        # the guard is on; absence of the key is the off switch (never a value).
        if 'coulomb_guard_w' in params:
            cgb = _to_numpy(params['coulomb_guard_b'])
            for i, a in enumerate(atoms):
                self._cg_b[a] = float(cgb[i])

        bk = _to_numpy(params['bond_k'])
        br0 = _to_numpy(params['bond_r0'])
        bis = [int(i) for i in data.bond_is]
        bjs = [int(j) for j in data.bond_js]
        for b in range(len(bis)):
            self._bonds[frozenset((atoms[bis[b]], atoms[bjs[b]]))] = (float(bk[b]), float(br0[b]))

        ak = _to_numpy(params['angle_k'])
        at0 = _to_numpy(params['angle_theta0'])
        ais = [int(i) for i in data.angle_is]
        ajs = [int(j) for j in data.angle_js]
        aks = [int(k) for k in data.angle_ks]
        for a in range(len(ais)):
            key = (atoms[ais[a]], atoms[ajs[a]], atoms[aks[a]])
            self._angles[key] = (float(ak[a]), float(at0[a]))

        pks = _to_numpy(params['proper_ks'])      # (n_propers, 6)
        pis = [int(i) for i in data.proper_is]
        pjs = [int(j) for j in data.proper_js]
        pks_i = [int(k) for k in data.proper_ks]
        pls = [int(l) for l in data.proper_ls]
        for p in range(len(pis)):
            key = (atoms[pis[p]], atoms[pjs[p]], atoms[pks_i[p]], atoms[pls[p]])
            self._propers[key] = [float(x) for x in pks[p]]

        iis = [int(i) for i in data.improper_is]
        ijs = [int(j) for j in data.improper_js]
        iks = [int(k) for k in data.improper_ks]
        ils = [int(l) for l in data.improper_ls]
        if 'oop_k' in params and len(iis):
            self.has_oop = True
            oop = _to_numpy(params['oop_k'])
            for m in range(len(iis)):
                key = (atoms[iis[m]], atoms[ijs[m]], atoms[iks[m]], atoms[ils[m]])
                self._impropers[key] = ('oop', float(oop[m]))
        elif len(iis):
            self.has_oop = False
            imk = _to_numpy(params['improper_ks'])   # (n_impropers, 2)
            for m in range(len(iis)):
                key = (atoms[iis[m]], atoms[ijs[m]], atoms[iks[m]], atoms[ils[m]])
                self._impropers[key] = ('periodic', [float(x) for x in imk[m]])

        gl = {
            'coulomb14scale': float(_to_numpy(params['coulomb14scale'])),
            'dexp_beta': float(_to_numpy(params['dexp_beta'])),
            'vdw14scale': float(_to_numpy(params.get('vdw14scale', 0.0))),
            'vdw13scale': float(_to_numpy(params.get('vdw13scale', 0.0))),
        }
        # The scalar repulsive-wall exponent exists only in the OLD (uniform-wall) form;
        # the per-atom-wall form drops it in favour of per-atom `dexp_b` (stored above).
        if 'dexp_alpha' in params:
            gl['dexp_alpha'] = float(_to_numpy(params['dexp_alpha']))
        # Coulomb-guard strength `w`, a scalar global, present iff the guard is on.
        if 'coulomb_guard_w' in params:
            gl['coulomb_guard_w'] = float(_to_numpy(params['coulomb_guard_w']))
        self._globals = gl

    # --------------------------------------------------------------- project
    def subset_for(self, all_atoms):
        '''
        Project the cached parameters onto a simulation's ``all_atoms`` array.

        ``all_atoms`` is the ordered :class:`chimerax.atomic.Atoms` from
        :class:`SimConstruct`; particle index = position in it. Terms with any
        atom outside ``all_atoms`` are dropped (the buffer-edge bonds into the
        non-simulated surroundings — always on frozen atoms, exactly as the
        AMBER path drops them via ``ignoreExternalBonds``).
        '''
        if not self.parameterised:
            raise RuntimeError('GarnetParameters.subset_for called before parameterise()')
        index = {a: i for i, a in enumerate(all_atoms)}
        n = len(all_atoms)
        charges = numpy.zeros(n)
        sigmas = numpy.zeros(n)
        epsilons = numpy.zeros(n)
        # Only materialise the optional per-atom arrays when this checkpoint predicts them,
        # so the builder can key off ``is None`` to choose the old vs new functional form.
        bee = numpy.zeros(n) if self._dexp_b else None
        cg_b = numpy.zeros(n) if self._cg_b else None
        missing = []
        for a, i in index.items():
            p = self._per_atom.get(a)
            if p is None:
                missing.append(a)
                continue
            charges[i], sigmas[i], epsilons[i] = p
            if bee is not None:
                bee[i] = self._dexp_b.get(a, 0.0)
            if cg_b is not None:
                cg_b[i] = self._cg_b.get(a, 0.0)

        bonds = []
        for key, (k, r0) in self._bonds.items():
            a, b = tuple(key)
            if a in index and b in index:
                bonds.append((index[a], index[b], k, r0))

        angles = []
        for (a, j, c), (k_theta, theta0) in self._angles.items():
            if a in index and j in index and c in index:
                angles.append((index[a], index[j], index[c], k_theta, theta0))

        propers = []
        for (a, b, c, d), ks in self._propers.items():
            if a in index and b in index and c in index and d in index:
                propers.append((index[a], index[b], index[c], index[d], ks))

        impropers = []
        for (a, b, c, d), (kind, payload) in self._impropers.items():
            if a in index and b in index and c in index and d in index:
                impropers.append((kind, (index[a], index[b], index[c], index[d]), payload))

        return SubsetParameters(n, charges, sigmas, epsilons, bonds, angles,
                                propers, impropers, self.has_oop, dict(self._globals), missing,
                                bee=bee, cg_b=cg_b)


import weakref
_PARAM_CACHE = weakref.WeakKeyDictionary()   # structure -> GarnetParameters


def get_garnet_parameters(structure, checkpoint_path=None, force=False, logger=None):
    '''
    Return a parameterised :class:`GarnetParameters` for ``structure``, caching
    it (keyed weakly to the model, so it is dropped when the model closes).

    Recomputes if the atom count changed since the last parameterisation (a
    proxy for a topology edit) or if ``force`` is True. Fragment-incremental
    reparameterisation is deferred to full Track B.
    '''
    cached = _PARAM_CACHE.get(structure)
    if (cached is not None and not force
            and cached.checkpoint_path == resolve_checkpoint_path(checkpoint_path)
            and len(cached._per_atom) == structure.num_atoms):
        return cached
    gp = GarnetParameters(structure, checkpoint_path=checkpoint_path)
    gp.parameterise(logger=logger)
    _PARAM_CACHE[structure] = gp
    return gp
