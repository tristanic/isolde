# @Author: Tristan Croll
# @Date:   02-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 02-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
MMFF94-based ligand parameterisation provider (pure-RDKit, no AmberTools).

The internal force terms (bonds, angles, torsions) come from MMFF94 via RDKit.
Partial charges are MMFF94 partial charges.  Per-atom LJ parameters are taken
from a GAFF2-derived element table so that Lorentz-Berthelot mixing with AMBER14
protein atoms is physically consistent.

This is a "best-effort" fallback for ligands that have no GARNET/espaloma/OpenFF
backend available.  Production simulations should prefer GARNET or an OpenFF-based
backend when accessible.
'''

import math
from .param_provider import LigandBackedParameterisationProvider

_KCAL_TO_KJ  = 4.184          # 1 kcal/mol = 4.184 kJ/mol
_ANG_TO_NM   = 0.1            # 1 Å = 0.1 nm
_DEG_TO_RAD  = math.pi / 180.0

# RDKit's MMFFMolProperties returns the *raw* MMFF94 force constants: bond kb in
# md/Å (millidyne/Å) and angle ka in md·Å/rad², NOT kcal/mol.  MMFF defines the
# harmonic energy as E = (kb/2)·143.9325·Δr² (and the analogous angle form), so
# this constant converts the raw force constant into kcal·mol⁻¹ per (Å² or rad²)
# before the usual kcal→kJ / Å→nm conversions.  Omitting it makes every bond and
# angle ~144x too weak.  (Torsion V1/V2/V3 are already in kcal/mol — no factor.)
_MDYNE_A_TO_KCAL = 143.9325

# Per-element LJ parameters (sigma/nm, epsilon/kJ·mol⁻¹) derived from GAFF2.
# Used for protein–ligand cross-term mixing via Lorentz-Berthelot rules.
_ELEMENT_LJ = {
    1:  (0.106, 0.0657),   # H
    6:  (0.339, 0.3598),   # C
    7:  (0.325, 0.7113),   # N
    8:  (0.296, 0.8786),   # O
    9:  (0.312, 0.2552),   # F
    15: (0.374, 0.8368),   # P
    16: (0.356, 1.0460),   # S
    17: (0.347, 1.1087),   # Cl
    35: (0.378, 1.6883),   # Br
    53: (0.418, 2.0920),   # I
}
_DEFAULT_LJ = (0.340, 0.50)  # fallback for uncommon elements


def _cx_residue_to_rdmol(residue):
    '''
    Build an RDKit RWMol from a ChimeraX residue and infer bond orders from
    3D geometry.  Requires explicit hydrogen atoms for reliable bond-order
    assignment.  ``charge=0`` is assumed.
    '''
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDetermineBonds

    cx_atoms = list(residue.atoms)
    atom_index = {a: i for i, a in enumerate(cx_atoms)}
    n = len(cx_atoms)

    mol = Chem.RWMol()
    conf = Chem.Conformer(n)
    for i, a in enumerate(cx_atoms):
        mol.AddAtom(Chem.Atom(a.element.number))
        c = a.coord
        conf.SetAtomPosition(i, (float(c[0]), float(c[1]), float(c[2])))
    mol.AddConformer(conf, assignId=True)

    bonds_obj = residue.atoms.intra_bonds
    a1_col, a2_col = bonds_obj.atoms
    for a1, a2 in zip(a1_col, a2_col):
        mol.AddBond(atom_index[a1], atom_index[a2], Chem.BondType.SINGLE)

    rdDetermineBonds.DetermineBondOrders(mol, charge=0)
    mol_ro = mol.GetMol()
    AllChem.AssignStereochemistryFrom3D(mol_ro)
    return mol_ro


def _mmff_system_to_ffxml(residue, mol) -> str:
    '''
    Extract MMFF94 parameters from *mol* and serialise them as a ForceField XML
    template compatible with ISOLDE's ForceField (USER_ prefix on registration).

    Force terms:
      - HarmonicBondForce : MMFF94 bond stretch
      - HarmonicAngleForce : MMFF94 angle bend
      - PeriodicTorsionForce : MMFF94 proper torsions (v1/v2/v3 decomposition)
      - NonbondedForce : MMFF94 partial charges + GAFF2-element LJ
    '''
    import xml.etree.ElementTree as ET
    from rdkit.Chem import rdForceFieldHelpers

    cx_atoms = list(residue.atoms)
    n_atoms  = mol.GetNumAtoms()
    rname    = residue.name
    type_names = [f'MMF_{rname}_{i}' for i in range(n_atoms)]

    # ---- MMFF94 properties (charges, atom types) ----------------------------
    props = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94')
    if props is None:
        raise RuntimeError(f'MMFF94 could not assign properties for residue {rname!r}')

    charges = [props.GetMMFFPartialCharge(i) for i in range(n_atoms)]

    def _lj(i):
        return _ELEMENT_LJ.get(mol.GetAtomWithIdx(i).GetAtomicNum(), _DEFAULT_LJ)

    sigma   = [_lj(i)[0] for i in range(n_atoms)]
    epsilon = [_lj(i)[1] for i in range(n_atoms)]

    root = ET.Element('ForceField')

    # AtomTypes
    at_el = ET.SubElement(root, 'AtomTypes')
    for i, cx_a in enumerate(cx_atoms):
        ET.SubElement(at_el, 'Type', {
            'name':    type_names[i],
            'class':   type_names[i],
            'element': cx_a.element.name.capitalize(),
            'mass':    str(cx_a.element.mass),
        })

    # Residue template
    res_el = ET.SubElement(ET.SubElement(root, 'Residues'), 'Residue', name=rname)
    for i, cx_a in enumerate(cx_atoms):
        ET.SubElement(res_el, 'Atom',
            name=cx_a.name, type=type_names[i], charge=str(charges[i]))
    bonds_obj = residue.atoms.intra_bonds
    a1_col, a2_col = bonds_obj.atoms
    for a1_name, a2_name in zip(a1_col.names, a2_col.names):
        ET.SubElement(res_el, 'Bond', atomName1=a1_name, atomName2=a2_name)

    # NonbondedForce (AMBER14-compatible 14-scaling)
    nb_el = ET.SubElement(root, 'NonbondedForce',
        coulomb14scale='0.8333333333333333', lj14scale='0.5')
    ET.SubElement(nb_el, 'UseAttributeFromResidue', name='charge')
    for i in range(n_atoms):
        ET.SubElement(nb_el, 'Atom',
            type=type_names[i], sigma=str(sigma[i]), epsilon=str(epsilon[i]))

    # ---- HarmonicBondForce --------------------------------------------------
    # GetMMFFBondStretchParams -> (bondType, kb, r0); kb in md/Å, r0 in Å.
    # MMFF94: E = (143.9325*kb/2)*Δr²  (kcal/mol, Δr in Å)
    # OpenMM: E = (k/2)*Δr²  with k in kJ/mol/nm²
    bf_el = None
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bp = props.GetMMFFBondStretchParams(mol, i, j)
        if bp is None:
            continue
        kb, r0 = bp[1], bp[2]
        k_kJ_nm2 = kb * _MDYNE_A_TO_KCAL * _KCAL_TO_KJ / (_ANG_TO_NM ** 2)
        r0_nm    = r0 * _ANG_TO_NM
        if bf_el is None:
            bf_el = ET.SubElement(root, 'HarmonicBondForce')
        ET.SubElement(bf_el, 'Bond',
            type1=type_names[i], type2=type_names[j],
            length=str(r0_nm), k=str(k_kJ_nm2))

    # ---- HarmonicAngleForce -------------------------------------------------
    # GetMMFFAngleBendParams -> (angleType, ka, theta0); ka in md·Å/rad²,
    # theta0 in degrees.  MMFF's harmonic term reduces to
    # E = (143.9325*ka/2)*Δθ_rad² (kcal/mol), so the same md->kcal constant
    # applies.  OpenMM: E = (k/2)*Δθ² with k in kJ/mol/rad², θ in rad.
    af_el = None
    for j in range(n_atoms):
        nbrs = [b.GetOtherAtomIdx(j)
                for b in mol.GetAtomWithIdx(j).GetBonds()]
        for ia in range(len(nbrs)):
            i = nbrs[ia]
            for ib in range(ia + 1, len(nbrs)):
                k = nbrs[ib]
                ap = props.GetMMFFAngleBendParams(mol, i, j, k)
                if ap is None:
                    continue
                ka, theta0 = ap[1], ap[2]
                ka_kJ   = ka * _MDYNE_A_TO_KCAL * _KCAL_TO_KJ
                th0_rad = theta0 * _DEG_TO_RAD
                if af_el is None:
                    af_el = ET.SubElement(root, 'HarmonicAngleForce')
                ET.SubElement(af_el, 'Angle',
                    type1=type_names[i], type2=type_names[j], type3=type_names[k],
                    angle=str(th0_rad), k=str(ka_kJ))

    # ---- PeriodicTorsionForce -----------------------------------------------
    # MMFF94: E = (V1/2)*(1+cos φ) + (V2/2)*(1-cos 2φ) + (V3/2)*(1+cos 3φ)
    # Mapped to OpenMM k*(1+cos(n*φ - phase)):
    #   V1 → k=V1/2, n=1, phase=0
    #   V2 → k=V2/2, n=2, phase=π
    #   V3 → k=V3/2, n=3, phase=0
    tf_el = None
    seen_torsions = set()
    for bond_jk in mol.GetBonds():
        j = bond_jk.GetBeginAtomIdx()
        k = bond_jk.GetEndAtomIdx()
        i_list = [b.GetOtherAtomIdx(j) for b in mol.GetAtomWithIdx(j).GetBonds()
                  if b.GetOtherAtomIdx(j) != k]
        l_list = [b.GetOtherAtomIdx(k) for b in mol.GetAtomWithIdx(k).GetBonds()
                  if b.GetOtherAtomIdx(k) != j]
        for i in i_list:
            for l in l_list:
                if l == i:
                    continue
                key = tuple(sorted([i, j, k, l]) + [j, k])
                if key in seen_torsions:
                    continue
                seen_torsions.add(key)
                tp = props.GetMMFFTorsionParams(mol, i, j, k, l)
                if tp is None:
                    continue
                # GetMMFFTorsionParams -> (torsionType, V1, V2, V3); V in kcal/mol.
                v1, v2, v3 = tp[1], tp[2], tp[3]
                if tf_el is None:
                    tf_el = ET.SubElement(root, 'PeriodicTorsionForce')
                for v, n_per, base_phase in ((v1, 1, 0.0),
                                              (v2, 2, math.pi),
                                              (v3, 3, 0.0)):
                    if abs(v) < 1e-6:
                        continue
                    k_kJ   = abs(v) * _KCAL_TO_KJ / 2.0
                    phase  = base_phase + (math.pi if v < 0.0 else 0.0)
                    ET.SubElement(tf_el, 'Proper',
                        type1=type_names[i], type2=type_names[j],
                        type3=type_names[k], type4=type_names[l],
                        periodicity1=str(n_per),
                        phase1=str(phase),
                        k1=str(k_kJ))

    return ET.tostring(root, encoding='unicode')


def _mmff_parameterise_residue(residue, logger) -> str:
    '''Run MMFF94 on the residue and return a ForceField XML string.'''
    from rdkit.Chem import AllChem

    mol = _cx_residue_to_rdmol(residue)
    result = AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    if result == -1:
        raise RuntimeError(
            f'MMFF94 could not set up a force field for {residue.name!r}. '
            f'The molecule may contain elements or bonding patterns not '
            f'supported by MMFF94.')
    return _mmff_system_to_ffxml(residue, mol)


class MMFFParameterisationProvider(LigandBackedParameterisationProvider):
    '''
    Parameterisation provider using AMBER14 for standard residues and MMFF94
    (via RDKit) for ligands.

    Requires only ``rdkit`` — no AmberTools, no OpenFF toolkit, no conda.

    Internal force constants (bonds, angles, torsions) come directly from
    MMFF94.  Per-atom LJ parameters are from a GAFF2-element table so that
    Lorentz-Berthelot mixing with AMBER14 protein atoms is physically
    consistent.  This is an approximation: for production work, prefer
    ``amber14+garnet`` or an OpenFF-based backend.
    '''

    @property
    def backend_forcefield_key(self):
        return 'amber14+mmff'

    def _parameterise_one(self, residue, logger):
        xml = _mmff_parameterise_residue(residue, logger)
        return xml, {}
