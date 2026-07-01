# @Author: Tristan Croll
# @Date:   01-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 01-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
GARNET-backed parameterisation provider.

:class:`GarnetParameterisationProvider` extends
:class:`ForceFieldParameterisationProvider`.  When ``sim_params.forcefield``
is ``'amber14+garnet'`` it intercepts :class:`UnparameterisedResiduesError`,
runs GARNET inference on the unmatched ligand residues, registers the resulting
XML templates under the ``USER_`` prefix, and retries the system build.

For any other forcefield value the class delegates transparently to the parent,
so installing this provider changes nothing for the pure AMBER14/CHARMM36 paths.

GARNET (Greener Group, MRC LMB, 2026) is a graph-neural-network forcefield
trained on quantum mechanical, condensed-phase and NMR data.  The model weights
are bundled with the ``garnetff`` package; no separate download is required.

Paper: arXiv:2603.16770 — "Training a force field for proteins and small
molecules from scratch" (Blanco-Gonzalez et al., 2026)

Dependencies (install once inside ChimeraX):
    pip install garnetff torch torch_geometric igraph

No openff.toolkit or RDKit is required — a self-contained shim topology is
built directly from ChimeraX atom/bond data.

VdW note
--------
GARNET uses a double-exponential potential internally.  This provider extracts
GARNET's per-atom sigma/epsilon values and registers them with a standard
Lennard-Jones ``NonbondedForce`` so that protein-ligand cross-term mixing
(Lorentz-Berthelot) works correctly through AMBER14.  The functional form
differs from GARNET's native potential, but the approximation is adequate for
interactive model building.
'''

from .param_provider import ForceFieldParameterisationProvider

# Atomic numbers supported by GARNET and their internal index/name (from garnetff source).
_ELEMENT_INDICES = {
    1: 0, 3: 1,  5: 2,  6: 3,  7: 4,  8: 5,  9: 6,
    11: 7, 12: 8, 14: 9, 15: 10, 16: 11,
    17: 12, 19: 13, 20: 14, 35: 15, 53: 16,
}
_ELEMENT_NAMES = [
    'H', 'Li', 'B', 'C', 'N', 'O', 'F',
    'Na', 'Mg', 'Si', 'P', 'S', 'Cl', 'K', 'Ca', 'Br', 'I',
]


# ---------------------------------------------------------------------------
# Shim topology classes
# Expose the interface expected by garnetff without requiring openff.toolkit.
# ---------------------------------------------------------------------------

class _GarnetFormalCharge:
    '''
    Minimal formal-charge quantity-like object.

    ``fc / fc.units`` and ``fc / fc.u`` both return the integer charge value,
    matching the two call sites inside garnetff source:
      - ``topology_to_data``:       ``int(a.formal_charge / a.formal_charge.units)``
      - ``find_identical_atoms``:   ``fci / fci.u``
    '''
    def __init__(self, value: int):
        self._v = int(value)
        self.units = 1
        self.u = 1

    def __truediv__(self, other):
        return self._v


class _AtomRef:
    '''Minimal atom reference used in angle/torsion index tuples.'''
    __slots__ = ('molecule_atom_index',)

    def __init__(self, idx: int):
        self.molecule_atom_index = idx


def _infer_formal_charge(cx_atom) -> int:
    '''Best-effort formal charge from IDATM type; 0 for the vast majority of atoms.'''
    t = cx_atom.idatm_type
    if t in ('N3+', 'N2+', 'O1+', 'Oar+'):
        return 1
    if t in ('O2-', 'S3-'):
        return -1
    return 0


def _garnet_element_name(atomic_number: int) -> str:
    idx = _ELEMENT_INDICES.get(atomic_number)
    if idx is None:
        raise ValueError(
            f'Atomic number {atomic_number} is not supported by GARNET. '
            'Supported elements: H Li B C N O F Na Mg Si P S Cl K Ca Br I.')
    return _ELEMENT_NAMES[idx]


def _compute_garnet_atom_names(cx_atoms) -> list:
    '''
    Generate GARNET-style element-count atom names in the same order as
    ``garnetff.data_to_xml`` would produce (C1, N1, O1, C2, …).
    '''
    counts = {}
    names = []
    for a in cx_atoms:
        el = _garnet_element_name(a.element.number)
        c = counts.get(el, 0) + 1
        counts[el] = c
        names.append(f'{el}{c}')
    return names


class _GarnetAtom:
    '''Wraps a ChimeraX atom to satisfy garnetff's atom interface.'''

    def __init__(self, cx_atom, idx: int):
        self._cx = cx_atom
        self._idx = idx
        self._fc = _GarnetFormalCharge(_infer_formal_charge(cx_atom))
        # Mutable dict — topology_to_openmm_system may overwrite residue_name for polymers.
        self.metadata = {
            'residue_name':   cx_atom.residue.name,
            'residue_number': cx_atom.residue.number,
            'insertion_code': cx_atom.residue.insertion_code or ' ',
        }

    @property
    def atomic_number(self) -> int:
        return self._cx.element.number

    @property
    def formal_charge(self) -> _GarnetFormalCharge:
        return self._fc

    @property
    def is_aromatic(self) -> bool:
        return self._cx.idatm_type in ('Car', 'Nar', 'Oar', 'Oar+', 'Sar')

    @property
    def bonded_atoms(self):
        '''Bonded atoms within the same residue (intra-residue neighbours).'''
        r = self._cx.residue
        return [nb for nb in self._cx.neighbors if nb.residue is r]


class _GarnetBond:
    __slots__ = ('atom1_index', 'atom2_index')

    def __init__(self, i: int, j: int):
        self.atom1_index = i
        self.atom2_index = j


class _GarnetMolecule:
    '''
    Shim molecule built from a single ChimeraX residue.

    Provides: n_atoms, atoms, bonds, angles, propers, amber_impropers,
    to_networkx(), to_smiles(), is_isomorphic_with().
    '''

    def __init__(self, residue):
        self._residue = residue
        cx_atoms = list(residue.atoms)
        self._cx_atoms = cx_atoms
        self._atoms = [_GarnetAtom(a, i) for i, a in enumerate(cx_atoms)]
        _id = {id(a): i for i, a in enumerate(cx_atoms)}

        # --- bonds (intra-residue) ---
        bonds = []
        seen_b = set()
        for a in cx_atoms:
            for nb in a.neighbors:
                if nb.residue is residue:
                    i, j = _id[id(a)], _id[id(nb)]
                    key = (min(i, j), max(i, j))
                    if key not in seen_b:
                        seen_b.add(key)
                        bonds.append(_GarnetBond(i, j))
        self._bonds = bonds

        # --- angles (all i–j–k paths through each central atom j) ---
        angles = []
        seen_a = set()
        for j_idx, j_cx in enumerate(cx_atoms):
            nbrs = [_id[id(nb)] for nb in j_cx.neighbors if nb.residue is residue]
            for ii in range(len(nbrs)):
                for ki in range(ii + 1, len(nbrs)):
                    i_idx, k_idx = nbrs[ii], nbrs[ki]
                    key = (min(i_idx, k_idx), j_idx, max(i_idx, k_idx))
                    if key not in seen_a:
                        seen_a.add(key)
                        angles.append([_AtomRef(i_idx), _AtomRef(j_idx), _AtomRef(k_idx)])
        self._angles = angles

        # --- propers (i–j–k–l paths across each bond j–k) ---
        propers = []
        seen_p = set()
        for b in bonds:
            j_idx, k_idx = b.atom1_index, b.atom2_index
            j_nbrs = [_id[id(nb)] for nb in cx_atoms[j_idx].neighbors
                      if nb.residue is residue and _id[id(nb)] != k_idx]
            k_nbrs = [_id[id(nb)] for nb in cx_atoms[k_idx].neighbors
                      if nb.residue is residue and _id[id(nb)] != j_idx]
            for i_idx in j_nbrs:
                for l_idx in k_nbrs:
                    if i_idx == l_idx:
                        continue
                    raw = (i_idx, j_idx, k_idx, l_idx)
                    key = min(raw, raw[::-1])
                    if key not in seen_p:
                        seen_p.add(key)
                        propers.append([_AtomRef(i_idx), _AtomRef(j_idx),
                                        _AtomRef(k_idx), _AtomRef(l_idx)])
        self._propers = propers

        # --- AMBER impropers (all combinations of 3 neighbours for each atom) ---
        from itertools import combinations
        impropers = []
        for j_idx, j_cx in enumerate(cx_atoms):
            nbrs = sorted(_id[id(nb)] for nb in j_cx.neighbors if nb.residue is residue)
            for i_idx, k_idx, l_idx in combinations(nbrs, 3):
                impropers.append([_AtomRef(j_idx), _AtomRef(i_idx),
                                   _AtomRef(k_idx), _AtomRef(l_idx)])
        self._impropers = impropers

    # ---- GARNET interface ----

    @property
    def n_atoms(self) -> int:
        return len(self._atoms)

    @property
    def atoms(self):
        return self._atoms

    @property
    def bonds(self):
        return self._bonds

    @property
    def angles(self):
        return self._angles

    @property
    def propers(self):
        return self._propers

    @property
    def amber_impropers(self):
        return self._impropers

    def to_networkx(self):
        import networkx as nx
        g = nx.Graph()
        for i, a in enumerate(self._atoms):
            g.add_node(i, atomic_number=a.atomic_number, formal_charge=a.formal_charge)
        for b in self._bonds:
            g.add_edge(b.atom1_index, b.atom2_index)
        return g

    def to_smiles(self) -> str:
        # Only used by garnetff to check for "." (multi-component molecule guard).
        # Residue names never contain ".".
        return self._residue.name

    def is_isomorphic_with(self, other) -> bool:
        # With a single-molecule topology this is never called (loop range is empty).
        return self is other


class _GarnetTopology:
    '''Single-molecule topology shim backed by a ChimeraX residue.'''

    def __init__(self, residue):
        self._mol = _GarnetMolecule(residue)

    @property
    def n_atoms(self) -> int:
        return self._mol.n_atoms

    @property
    def n_molecules(self) -> int:
        return 1

    @property
    def molecules(self):
        return [self._mol]

    @property
    def box_vectors(self):
        return None

    def to_openmm(self):
        '''
        Build an OpenMM Topology with GARNET-style element-count atom names
        (C1, N1, O1, C2, …) — these must match the names in the GARNET-
        generated ForceField XML that ``createSystem`` will be called with.
        '''
        from openmm.app import Topology as OMMTopology, Element as OMMElement
        top = OMMTopology()
        chain = top.addChain()
        res = top.addResidue(self._mol._residue.name, chain)
        garnet_names = _compute_garnet_atom_names(self._mol._cx_atoms)
        omm_atoms = []
        for ga, gname in zip(self._mol._atoms, garnet_names):
            elem = OMMElement.getByAtomicNumber(ga.atomic_number)
            omm_atoms.append(top.addAtom(gname, elem, res))
        for b in self._mol._bonds:
            top.addBond(omm_atoms[b.atom1_index], omm_atoms[b.atom2_index])
        return top


# ---------------------------------------------------------------------------
# System → FFXML conversion
# ---------------------------------------------------------------------------

def _garnet_system_to_ffxml(residue, system) -> str:
    '''
    Convert a GARNET-produced OpenMM System to a ForceField XML template
    string suitable for
    ``ForceField.loadFile(io.StringIO(xml), resname_prefix='USER_')``.

    Atom indices in *system* correspond to ``list(residue.atoms)`` order
    (guaranteed by :class:`_GarnetTopology` construction).

    Charge/LJ strategy
    ------------------
    GARNET stores partial charges via ``<UseAttributeFromResidue name="charge"/>``
    in its NonbondedForce and uses a CustomNonbondedForce for VdW.  After
    ``createSystem()`` the NonbondedForce carries per-particle charges with
    sigma=1/epsilon=0, and the CustomNonbondedForce carries the real sigma/epsilon.

    This function reads charges from NonbondedForce and sigma/epsilon from
    CustomNonbondedForce, then writes a single AMBER14-compatible NonbondedForce
    entry combining all three — allowing correct Lorentz-Berthelot mixing with
    protein atoms parameterised by AMBER14.
    '''
    import xml.etree.ElementTree as ET
    from openmm import unit as omm_unit
    from openmm.openmm import (HarmonicBondForce, HarmonicAngleForce,
                               PeriodicTorsionForce, NonbondedForce,
                               CustomNonbondedForce)

    cx_atoms = list(residue.atoms)
    n_atoms = system.getNumParticles()
    rname = residue.name
    type_names = [f'GRN_{rname}_{i}' for i in range(n_atoms)]

    # --- extract non-bonded parameters ---
    charges = [0.0] * n_atoms
    sigma = [0.35] * n_atoms    # fallback: generic sigma in nm
    epsilon = [0.40] * n_atoms  # fallback: generic epsilon in kJ/mol

    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            for ni in range(force.getNumParticles()):
                q, _sig, _eps = force.getParticleParameters(ni)
                charges[ni] = q.value_in_unit(omm_unit.elementary_charge)
        elif isinstance(force, CustomNonbondedForce):
            for ni in range(force.getNumParticles()):
                params = force.getParticleParameters(ni)
                # GARNET's CustomNonbondedForce defines PerParticleParameters
                # "sigma" and "epsilon" in that order (nm, kJ/mol).
                if len(params) >= 2:
                    sigma[ni] = float(params[0])
                    epsilon[ni] = float(params[1])

    root = ET.Element('ForceField')

    # --- AtomTypes ---
    at_el = ET.SubElement(root, 'AtomTypes')
    for i, cx_a in enumerate(cx_atoms):
        ET.SubElement(at_el, 'Type',
            name=type_names[i],
            cls=type_names[i],
            element=cx_a.element.name.capitalize(),
            mass=str(cx_a.element.mass))

    # --- Residue template ---
    res_el = ET.SubElement(ET.SubElement(root, 'Residues'), 'Residue', name=rname)
    for i, cx_a in enumerate(cx_atoms):
        ET.SubElement(res_el, 'Atom',
            name=cx_a.name,
            type=type_names[i],
            charge=str(charges[i]))
    # Bond connectivity (atom names, not parameters)
    bonds = residue.atoms.intra_bonds
    a1_col, a2_col = bonds.atoms
    for a1_name, a2_name in zip(a1_col.names, a2_col.names):
        ET.SubElement(res_el, 'Bond', atomName1=a1_name, atomName2=a2_name)

    # --- NonbondedForce (AMBER14-compatible 14-scaling) ---
    nb_el = ET.SubElement(root, 'NonbondedForce',
        coulomb14scale='0.8333333333333333', lj14scale='0.5')
    ET.SubElement(nb_el, 'UseAttributeFromResidue', name='charge')
    for i in range(n_atoms):
        ET.SubElement(nb_el, 'Atom',
            type=type_names[i],
            sigma=str(sigma[i]),
            epsilon=str(epsilon[i]))

    # --- Bonded forces ---
    for force in system.getForces():
        if isinstance(force, HarmonicBondForce):
            bf_el = ET.SubElement(root, 'HarmonicBondForce')
            for bi in range(force.getNumBonds()):
                i1, i2, r0, k = force.getBondParameters(bi)
                ET.SubElement(bf_el, 'Bond',
                    type1=type_names[i1], type2=type_names[i2],
                    length=str(r0.value_in_unit(omm_unit.nanometer)),
                    k=str(k.value_in_unit(
                        omm_unit.kilojoule_per_mole / omm_unit.nanometer**2)))

        elif isinstance(force, HarmonicAngleForce):
            af_el = ET.SubElement(root, 'HarmonicAngleForce')
            for ai in range(force.getNumAngles()):
                i1, i2, i3, theta0, k = force.getAngleParameters(ai)
                ET.SubElement(af_el, 'Angle',
                    type1=type_names[i1], type2=type_names[i2],
                    type3=type_names[i3],
                    angle=str(theta0.value_in_unit(omm_unit.radian)),
                    k=str(k.value_in_unit(
                        omm_unit.kilojoule_per_mole / omm_unit.radian**2)))

        elif isinstance(force, PeriodicTorsionForce):
            tf_el = ET.SubElement(root, 'PeriodicTorsionForce')
            for ti in range(force.getNumTorsions()):
                i1, i2, i3, i4, per, phase, k = force.getTorsionParameters(ti)
                # All torsion entries are written as <Proper>.  Improper out-of-
                # plane terms from GARNET cannot be distinguished at system level;
                # those that have no bonded 1-2-3-4 path will be silently skipped
                # by OpenMM's template matching, which is acceptable for model
                # building.
                ET.SubElement(tf_el, 'Proper',
                    type1=type_names[i1], type2=type_names[i2],
                    type3=type_names[i3], type4=type_names[i4],
                    periodicity1=str(per),
                    phase1=str(phase.value_in_unit(omm_unit.radian)),
                    k1=str(k.value_in_unit(omm_unit.kilojoule_per_mole)))

    return ET.tostring(root, encoding='unicode')


# ---------------------------------------------------------------------------
# Provider
# ---------------------------------------------------------------------------

class GarnetParameterisationProvider(ForceFieldParameterisationProvider):
    '''
    Parameterisation provider that uses AMBER14 for standard residues and
    GARNET for ligands that have no existing template.

    When ``sim_params.forcefield`` is ``'amber14'`` or ``'charmm36'`` the class
    behaves identically to :class:`ForceFieldParameterisationProvider`.
    When it is ``'amber14+garnet'`` it intercepts
    :class:`UnparameterisedResiduesError`, runs GARNET on each unmatched
    non-polymer residue, registers the XML template under the ``USER_`` prefix,
    and retries.

    Parameters
    ----------
    forcefield_mgr : ForcefieldMgr
    '''

    def __init__(self, forcefield_mgr):
        super().__init__(forcefield_mgr)
        # Session-lifetime cache: {residue_name: xml_str}
        self._template_cache = {}

    def build_system(self, sim_construct, sim_params, logger=None):
        if sim_params.forcefield != 'amber14+garnet':
            return super().build_system(sim_construct, sim_params, logger)
        ff  = self._forcefield_mgr['amber14']
        ldb = self._forcefield_mgr.ligand_db('amber14')
        return self._build_with_garnet_fallback(
            ff,
            sim_construct.all_atoms,
            sim_construct.all_residues,
            sim_params,
            ldb,
            logger)

    def _build_with_garnet_fallback(self, ff, atoms, all_residues,
                                    sim_params, ligand_db, logger):
        '''
        3-attempt retry loop:
          0 → stale-override recovery (same as parent)
          1 → GARNET fill-in for genuinely unparameterised ligands
          2 → propagate error
        '''
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
            return top, system

    def _parameterise_unmatched(self, unmatched_residues, forcefield, logger):
        '''
        Run GARNET on each unmatched non-polymer residue and register the
        resulting template under the USER_ prefix.
        '''
        import io
        from chimerax.atomic import Residue
        for r in unmatched_residues:
            if r.polymer_type != Residue.PT_NONE:
                continue
            name = r.name
            if name in self._template_cache:
                xml = self._template_cache[name]
            else:
                if logger is not None:
                    logger.status(f'GARNET: parameterising {name}...')
                xml = self._garnet_parameterise_one(r, logger)
                self._template_cache[name] = xml
                if logger is not None:
                    logger.info(f'GARNET parameterisation complete for {name}.')
            forcefield.loadFile(io.StringIO(xml), resname_prefix='USER_')

    def _garnet_parameterise_one(self, residue, logger):
        '''
        Build shim topology → GARNET inference → ISOLDE-compatible XML string.

        The shim topology exposes the interface garnetff expects without
        requiring openff.toolkit or RDKit.
        '''
        from garnetff import garnet
        from openmm.app import NoCutoff

        topology = _GarnetTopology(residue)
        system, _top_omm = garnet.topology_to_openmm_system(
            topology,
            nonbondedMethod=NoCutoff,
            constraints=None)

        return _garnet_system_to_ffxml(residue, system)
