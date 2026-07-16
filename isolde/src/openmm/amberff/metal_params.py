# @Author: Tristan Croll
# @Date:   14-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 14-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Empirical metal-site parameters for ISOLDE's metal-coordination MD
parameterisation (the metal analogue of ``boundary_params``).

Metal-containing ligands and coordination sites (hemes, chlorophylls, zinc
fingers, ...) cannot be typed by ANTECHAMBER/GAFF2 -- GAFF has no transition-metal
atom types. The reference AMBER workflow, MCPB.py, builds a *bonded* metal model
whose metal-ligand force constants come from a QM Hessian (Seminario method) and
whose charges come from RESP -- requiring Gaussian/GAMESS/ORCA and hours of expert
effort.

ISOLDE targets a different regime: hold a *sensible coordination geometry* while
the user rebuilds into density. That does not need QM. We produce the same *kind*
of artifact MCPB does (a bonded metal-site template: metal-ligand bonds +
coordination angles + redistributed charges), but:

* equilibrium **distances** come from a curated table (seeded from crystallographic
  averages and the shipped bonded templates ``BRYCE_HEM``/``iron_sulfur.xml``), with
  a UFF-radii fallback for the long tail;
* equilibrium **angles** are obtained by snapping each *observed* donor-metal-donor
  angle to the nearest canonical value (90/109.5/120/180 deg) -- which adapts to
  octahedral / tetrahedral / square-planar geometry straight from the model without
  a rigid polyhedron/oxidation-state table;
* **force constants** are soft, in the range of the shipped bonded templates (bond
  ~40-80k kJ/mol/nm^2, angle ~200-420 kJ/mol/rad^2) -- stiff enough to hold
  coordination, soft enough that a mis-placed donor relaxes rather than exploding.
  This is deliberately *not* QM-accurate; it is fit for interactive model building,
  not free-energy work.
* the metal's formal (oxidation-state) charge is spread partly onto the donor atoms
  (:func:`metal_charge_split`), a cheap approximation of the charge transfer a
  QM/RESP fit would capture (real metal charges sit well below the oxidation state).

All quantities here are in **OpenMM ffXML units**: lengths nm, angles radians,
bond k kJ/mol/nm^2, angle k kJ/mol/rad^2. These functions are pure (no ChimeraX /
OpenMM imports) so they are cheap to unit-test headlessly.
'''

import math

#: Ultimate fallback bond force constant (kJ/mol/nm^2), used only when a pair's
#: reduced mass is unknown so the frequency estimator below cannot run.
DEFAULT_BOND_K = 62760.0        # ~150 kcal/mol/A^2

#: Coordination ANGLE force constant (kJ/mol/rad^2). Set at MCPB's validation ceiling
#: for metal-ligand angles (100 kcal/mol/rad^2 ~= 418 kJ/mol/rad^2) and the top of the
#: shipped ``BRYCE_HEM`` range. The angles -- not the (physically soft) dative bonds --
#: are the lever that HOLDS a lone axial donor (e.g. a heme His) on the coordination
#: axis during interactive fitting: the bond fixes the metal-donor distance, the angle
#: terms fix the donor's *direction*. A physically faithful (weak) axial bond would let
#: the donor swing at fixed distance unless the angles are firm, so we deliberately
#: keep the bonds at their measured stiffness and hold geometry through stiffer angles.
DEFAULT_ANGLE_K = 418.4         # 100 kcal/mol/rad^2 (MCPB angle ceiling)

#: Physical constant for the harmonic-oscillator relation k = mu*(2*pi*c*nu)^2, in the
#: mixed units used here: k[kJ/mol/nm^2] = _K_FREQ_CONST * mu[amu] * nu[cm^-1]^2.
#: Derivation: (N/m -> kJ/mol/nm^2 = 602.214) * (amu -> kg = 1.66054e-27) *
#: (2*pi*c[cm/s])^2 with c = 2.99792e10.
_K_FREQ_CONST = 0.0354815

#: Fallback metal-donor stretch frequency (cm^-1) for pairs with no measured value.
DEFAULT_DONOR_FREQ = 300.0

#: Atomic masses (amu) for the reduced-mass term of the frequency estimator.
ATOMIC_MASS = {
    'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.06, 'Se': 78.97,
    'P': 30.974, 'F': 18.998, 'Cl': 35.45, 'Br': 79.904, 'I': 126.904,
    'Li': 6.94, 'Na': 22.990, 'K': 39.098, 'Rb': 85.468, 'Cs': 132.905,
    'Mg': 24.305, 'Ca': 40.078, 'Sr': 87.62, 'Ba': 137.327,
    'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546,
    'Zn': 65.38, 'Cd': 112.414, 'Hg': 200.592,
    'Mo': 95.95, 'W': 183.84, 'V': 50.942, 'Cr': 51.996, 'Ru': 101.07, 'Al': 26.982,
}

#: Multiplier applied to the GAFF planar-improper force constants on the atoms of a
#: metal-coordinating conjugated macrocycle (porphyrin/chlorin/corrin ring core).
#: GAFF gives aromatic out-of-plane impropers only ~1.1 kcal/mol (~4.6 kJ/mol);
#: over a 24-atom porphyrin those errors accumulate and the ring visibly puckers/
#: wobbles around the metal. Stiffening just the macrocycle impropers holds the near-
#: planar geometry (real hemes are flat bar a slight dome) without touching the
#: substituent impropers (carboxylate/vinyl), which GAFF already gets right. Applies
#: ONLY to the ring-system core selected by ring membership + metal-donor seeding, so
#: flexible pendants keep their normal flexibility. Tunable -- raise to flatten
#: further, lower toward 1.0 to restore stock GAFF behaviour.
MACROCYCLE_IMPROPER_SCALE = 4.0

#: Urey-Bradley 1-3 force constant (kJ/mol/nm^2) for a polynuclear inorganic cluster
#: core (iron-sulfur, ...). These stiff second-neighbour (metal...metal / bridge...bridge)
#: distance terms -- not the angles -- are what stop the cage collapsing; the shipped
#: hand-curated ``iron_sulfur.xml`` uses k=100000 with its equilibrium distances taken
#: from an idealised cluster. We reuse that k and take the equilibrium distance from the
#: MODEL (observed), so an arbitrary cluster is frozen at its own reliable geometry.
UREY_BRADLEY_K = 100000.0

#: Lennard-Jones (sigma nm, epsilon kJ/mol) for a bare bridging inorganic core atom
#: (mu-sulfide / mu-oxide), which never goes through GAFF/AM1-BCC. Sulfide seeded from
#: ``iron_sulfur.xml`` (its ST/SB type); oxide from an ff14SB carboxylate-oxygen-like
#: value. Keyed by element symbol.
CORE_ATOM_LJ = {
    'S': (0.4098134103445561, 1.046),      # iron_sulfur.xml ST/SB
    'O': (0.2959921901925, 0.87864),       # ff14SB O2-like
}

#: Idealised cluster geometry as reference atom COORDINATES (Angstrom, arbitrary frame),
#: taken from geostd's monomer entries -- which carry a real high-quality structure's
#: coordinates, not a computed idealisation. Keyed by CCD id -> {atom_name: (x, y, z)}.
#: The cluster-internal bond lengths, angles and Urey-Bradley 1-3 distances are computed
#: from these, NOT from the modelled coordinates (often badly distorted at low
#: resolution) and NOT from the CCD 'ideal' coordinates (geometrically unreliable for
#: metal clusters -- F3S comes out with ~83/147 deg S-Fe-S angles instead of ~104/113).
#: The model's cluster is mapped onto these coordinates by GRAPH ISOMORPHISM (element +
#: connectivity), NOT by atom name -- cluster atom labels (which sulfur is "S1") are
#: permuted between structures. Extend as more clusters are validated; an uncurated
#: cluster falls back to the highest-resolution PDB instance
#: (covalent._rcsb_cluster_coords), then to the modelled geometry.
CLUSTER_IDEAL_COORDS = {
    'F3S': {
        'FE1': (-36.7440, -10.7330, 10.6830), 'FE3': (-36.9320, -8.7940, 12.5680),
        'FE4': (-34.8430, -10.5140, 12.5870), 'S1': (-34.7520, -8.3090, 13.0650),
        'S2': (-37.0670, -11.0590, 12.9440), 'S3': (-34.4830, -11.0280, 10.4120),
        'S4': (-37.3730, -8.6100, 10.3930),
    },
    'SF4': {
        'FE1': (16.2190, -0.9470, 41.5710), 'FE2': (14.9120, 1.3570, 42.1980),
        'FE3': (17.2470, 0.7470, 43.4250), 'FE4': (14.9620, -0.7040, 43.9540),
        'S1': (15.3010, 1.5110, 44.4140), 'S2': (17.0770, -1.5320, 43.5500),
        'S3': (13.9710, -0.7190, 41.9620), 'S4': (16.9920, 1.1960, 41.2660),
    },
    'FES': {
        'FE1': (-0.0125, 0.2368, -1.3376), 'FE2': (0.0125, -0.0778, 1.3376),
        'S1': (1.6915, -0.2992, -0.0603), 'S2': (-1.6915, 0.4582, 0.0603),
    },
}

#: Curated metal-donor ``(equilibrium bond length nm, target stretch frequency
#: cm^-1)``, keyed by ``(metal, donor)`` -- or ``(metal, donor, 'axial')`` for a
#: weaker axial/dative bond where that distinction matters. Both element symbols are
#: capitalised as ChimeraX reports them (``'Fe'``, ``'Zn'``, ``'N'``).
#:
#: The force constant is NOT stored -- it is DERIVED from the frequency via the
#: harmonic relation ``k = mu*(2*pi*c*nu)^2`` (:func:`bond_k_from_frequency`), so each
#: bond's stiffness follows the *measured* metal-ligand stretch rather than a flat
#: guess. Sanity check: Fe-N(porphyrin) at 340 cm^-1 gives ~46000 kJ/mol/nm^2, which
#: matches the shipped ``BRYCE_HEM`` Fe-NP and the observed equatorial Fe-N mode.
#:
#: Frequencies grounded in resonance-Raman / far-IR data where noted; others are
#: representative. Distances are crystallographic averages. This is the place to
#: correct a stiffness that matters -- adjust the FREQUENCY (the physical quantity),
#: not an opaque k.
#:
#: Axial note: the axial His-Fe bond of a heme is physically weaker than the
#: equatorial pyrrole bonds -- ~300 cm^-1 for a 6-coordinate low-spin centre
#: (cytochromes; the common case), down to ~220 cm^-1 for 5-coordinate high-spin
#: (deoxy globins). We default to the low-spin value; a lone axial donor is HELD by
#: the (stiff) coordination angles, not by making this bond artificially strong.
METAL_SITE_PARAMS = {
    ('Fe', 'N'): (0.201, 340.0),          # heme pyrrole (equatorial); cf. BRYCE_HEM Fe-NP
    ('Fe', 'N', 'axial'): (0.201, 300.0), # axial His, low-spin (cyt c); high-spin ~220
    ('Fe', 'O'): (0.205, 400.0),
    ('Fe', 'S'): (0.230, 350.0),          # Fe-S(Cys), cf. rubredoxin 314-376 cm^-1
    ('Zn', 'N'): (0.205, 240.0),          # Zn-N(His) imidazole ~200-260 cm^-1
    ('Zn', 'O'): (0.203, 330.0),          # Asp/Glu/water
    ('Zn', 'S'): (0.231, 320.0),          # Cys thiolate
    ('Mg', 'N'): (0.220, 300.0),
    ('Mg', 'O'): (0.209, 350.0),
    ('Mn', 'N'): (0.221, 260.0),
    ('Mn', 'O'): (0.218, 320.0),
    ('Ca', 'O'): (0.240, 280.0),          # Ca-O soft / labile
    ('Cu', 'N'): (0.202, 290.0),
    ('Cu', 'S'): (0.230, 320.0),
    ('Ni', 'N'): (0.210, 300.0),
    ('Co', 'N'): (0.210, 300.0),
    ('Na', 'O'): (0.240, 230.0),
    ('K',  'O'): (0.280, 200.0),
}

#: UFF valence-bond radii (Angstrom) -- the long-tail fallback for metal-donor
#: equilibrium distance when the pair is absent from :data:`METAL_SITE_PARAMS`.
#: r0 ~= (r_metal + r_donor). Values from Rappe et al. 1992 (UFF), covering the
#: common donors plus a spread of metals. Approximate by design -- the curated
#: table above should carry any pair that matters.
UFF_RADII = {
    'H': 0.354, 'C': 0.757, 'N': 0.700, 'O': 0.658, 'S': 1.064, 'Se': 1.190,
    'F': 0.668, 'Cl': 1.044, 'Br': 1.192, 'I': 1.382, 'P': 1.117,
    'Li': 1.336, 'Na': 1.539, 'K': 1.953, 'Rb': 2.260, 'Cs': 2.570,
    'Mg': 1.421, 'Ca': 1.761, 'Sr': 1.914, 'Ba': 2.260,
    'Mn': 1.382, 'Fe': 1.412, 'Co': 1.241, 'Ni': 1.164, 'Cu': 1.302,
    'Zn': 1.193, 'Cd': 1.424, 'Hg': 1.340,
    'Mo': 1.484, 'W': 1.526, 'V': 1.402, 'Cr': 1.345, 'Ru': 1.478,
}

#: Canonical coordination angles (radians): linear-ish 90, tetrahedral 109.47,
#: trigonal 120, and trans/linear 180. An observed donor-metal-donor angle is
#: snapped to whichever of these it is closest to (see :func:`snap_angle`), so the
#: emitted geometry cleans up toward the ideal polyhedron the model implies without
#: needing an explicit coordination-number / oxidation-state classifier.
CANONICAL_ANGLES = (math.radians(90.0), math.radians(109.47),
                    math.radians(120.0), math.radians(180.0))

#: Fallback oxidation state per metal element when the modelled atom carries no
#: formal charge (used only to seed :func:`metal_charge_split`).
DEFAULT_OX_STATE = {
    'Li': 1, 'Na': 1, 'K': 1, 'Rb': 1, 'Cs': 1, 'Cu': 2, 'Ag': 1,
    'Mg': 2, 'Ca': 2, 'Sr': 2, 'Ba': 2, 'Zn': 2, 'Cd': 2, 'Hg': 2,
    'Mn': 2, 'Fe': 2, 'Co': 2, 'Ni': 2, 'Ca2': 2,
    'Fe2': 2, 'Fe3': 3,
    'Al': 3, 'Cr': 3,
}


#: ff14SB atom type for a metal-coordinating (deprotonated) donor atom of a standard
#: residue. This lets the metal-donor bond/angle parameters be keyed by atom TYPE
#: without emitting a custom template for the coordinating residue -- it keeps its
#: normal ISOLDE template. The coordinating N of a His is the DEProtonated ring N,
#: which is ff14SB type ``NB`` in both HID and HIE, so the protonation variant is
#: irrelevant. Keyed ``(residue_name, atom_name)`` with an element fallback.
_COORD_DONOR_TYPE = {
    ('HIS', 'ND1'): 'NB', ('HIS', 'NE2'): 'NB',
    ('HID', 'ND1'): 'NB', ('HID', 'NE2'): 'NB',
    ('HIE', 'ND1'): 'NB', ('HIE', 'NE2'): 'NB',
    ('HIP', 'ND1'): 'NB', ('HIP', 'NE2'): 'NB',
    ('CYS', 'SG'): 'SH', ('CYM', 'SG'): 'SH', ('CYX', 'SG'): 'S',
    ('MET', 'SD'): 'S',
    ('ASP', 'OD1'): 'O2', ('ASP', 'OD2'): 'O2',
    ('GLU', 'OE1'): 'O2', ('GLU', 'OE2'): 'O2',
    ('SER', 'OG'): 'OH', ('THR', 'OG1'): 'OH', ('TYR', 'OH'): 'OH',
    ('ASN', 'OD1'): 'O', ('GLN', 'OE1'): 'O',
}
_COORD_DONOR_ELEMENT_FALLBACK = {'N': 'NB', 'O': 'O', 'S': 'S'}


def coordinating_donor_type(residue_name, atom_name, element):
    '''ff14SB atom type for a standard-residue coordinating donor, or ``None`` if the
    element is not a recognised donor.'''
    t = _COORD_DONOR_TYPE.get((residue_name, atom_name))
    if t is not None:
        return t
    return _COORD_DONOR_ELEMENT_FALLBACK.get(element)


def bond_k_from_frequency(metal_element, donor_element, nu_cm1):
    '''Harmonic bond force constant (kJ/mol/nm^2) for a metal-donor stretch of
    frequency ``nu_cm1`` (cm^-1), from ``k = mu*(2*pi*c*nu)^2`` where ``mu`` is the
    two-atom reduced mass. This is the clean physical basis: the force constant
    follows the *measured* bond stiffness. Falls back to :data:`DEFAULT_BOND_K` if a
    mass is unknown (so the estimator can never crash a parameterisation).'''
    m1 = ATOMIC_MASS.get(metal_element)
    m2 = ATOMIC_MASS.get(donor_element)
    if m1 is None or m2 is None:
        return DEFAULT_BOND_K
    mu = m1 * m2 / (m1 + m2)
    return _K_FREQ_CONST * mu * nu_cm1 * nu_cm1


def metal_bond_params(metal_element, donor_element, axial=False):
    '''``(r0_nm, k)`` for a metal-donor bond. The distance comes from the curated
    table (crystallographic average) or the UFF-radii sum; the force constant is
    DERIVED from the pair's target stretch frequency via
    :func:`bond_k_from_frequency` (not a stored constant). Returns ``None`` only when
    even the UFF radii are unknown for one of the elements.

    ``axial=True`` selects the weaker axial/dative frequency where the table defines
    one (currently Fe-N: the heme His bond). It is used for donors contributed by a
    separate standard residue -- the heme-His case -- and is a no-op for pairs with no
    axial entry, so it never changes an equatorial / in-ligand bond.'''
    r0 = freq = None
    if axial:
        entry = METAL_SITE_PARAMS.get((metal_element, donor_element, 'axial'))
        if entry is not None:
            r0, freq = entry
    if r0 is None:
        entry = METAL_SITE_PARAMS.get((metal_element, donor_element))
        if entry is not None:
            r0, freq = entry
    if r0 is None:
        rm = UFF_RADII.get(metal_element)
        rd = UFF_RADII.get(donor_element)
        if rm is None or rd is None:
            return None
        r0, freq = (rm + rd) * 0.1, DEFAULT_DONOR_FREQ
    return (r0, bond_k_from_frequency(metal_element, donor_element, freq))


def snap_angle(observed_radians):
    '''Snap an observed donor-metal-donor angle to the nearest canonical
    coordination angle (see :data:`CANONICAL_ANGLES`).'''
    return min(CANONICAL_ANGLES, key=lambda a: abs(a - observed_radians))


def angle_k():
    '''Force constant (kJ/mol/rad^2) for a coordination angle -- firm by design (see
    :data:`DEFAULT_ANGLE_K`): the angles hold coordination geometry, especially a lone
    axial donor, so its physically-soft dative bond need not be over-stiffened.'''
    return DEFAULT_ANGLE_K


def metal_charge_split(ox_state, n_donors, keep_fraction=0.5):
    '''Approximate the charge transfer a QM/RESP fit would capture: keep a fraction
    of the metal's formal (oxidation-state) charge on the metal and spread the
    remainder evenly over its ``n_donors`` coordinating atoms.

    Real metal partial charges sit well below the oxidation state (MCPB's own
    validation notes the metal RESP charge is "often < +1"); ``keep_fraction=0.5``
    is a deliberately simple default. The split conserves total charge: the metal's
    contribution to the site is ``q_metal + n_donors * per_donor_delta ==
    ox_state``.

    Returns ``(q_metal, per_donor_delta)``.
    '''
    q_metal = keep_fraction * ox_state
    if n_donors <= 0:
        return ox_state, 0.0
    per_donor_delta = (ox_state - q_metal) / n_donors
    return q_metal, per_donor_delta


def guess_ox_state(element, modelled_formal_charge=None):
    '''Best guess at a metal's oxidation state: the modelled formal charge if it is
    non-zero, else the per-element default, else +2.'''
    if modelled_formal_charge:
        return int(modelled_formal_charge)
    return DEFAULT_OX_STATE.get(element, 2)
