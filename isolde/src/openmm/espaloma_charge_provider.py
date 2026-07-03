# @Author: Tristan Croll
# @Date:   03-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 03-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Hybrid ligand parameterisation: espaloma-charge for partial charges +
MMFF94 (RDKit) for bonds/angles/torsions + GAFF2-element LJ table.

This is a pip-installable alternative to full espaloma::

    pip install espaloma-charge dgl

No conda or AmberTools required.

- ``espaloma-charge`` (0.0.8+) provides GNN-based partial charges trained to
  reproduce AM1-BCC quality.
- All bonded parameters (bonds, angles, torsions) come from MMFF94 via RDKit.
- LJ parameters use the same GAFF2-element table as ``amber14+mmff`` so that
  Lorentz-Berthelot mixing with AMBER14 protein atoms is physically consistent.

DGL compatibility
-----------------
``espaloma-charge`` pulls in DGL, whose ``graphbolt`` sampling subsystem ships
version-matched native libraries that do not exist for recent PyTorch builds
(e.g. DGL 2.2.1 tops out at torch 2.3, but ChimeraX ships torch >= 2.12).  On
such installs a bare ``import dgl`` dies inside ``graphbolt`` before core DGL
loads.  :func:`_ensure_dgl_importable` neutralises that at import time so the
backend works out of the box without hand-patching the DGL install.

The model checkpoint (a pickled ``nn.Module``) is loaded once per session via
:func:`_load_espaloma_model` with ``weights_only=False`` (PyTorch >= 2.6 defaults
it to ``True`` and would otherwise refuse the checkpoint), and cached under
``torch.hub`` rather than the current working directory that stock
``espaloma_charge.charge`` uses.

Known limitations
-----------------
- **Charged ligands are not supported.**  Bond orders are inferred from 3D
  geometry via ``DetermineBondOrders(mol, charge=0)`` (in
  :func:`.mmff_provider._cx_residue_to_rdmol`), i.e. the ligand is assumed
  neutral.  A charged species such as GTP (net -8) fails bond perception with
  "Final molecular charge does not match input"; the residue then falls back to
  AMBER14.  This is a safe, expected fallback, not a crash.
- **Bonded terms are MMFF94, not a self-consistent AMBER/GAFF field.**  Charges
  are AM1-BCC-quality, but the bonded/VdW terms are an approximation; prefer
  ``amber14+garnet`` for charge-sensitive production work.
- **Good smoke-test ligands** (neutral, MMFF94-supported, unambiguous bonds):
  glycerol (``GOL``), benzene (``BNZ``), ethylene glycol (``EDO``).
'''

import math
from .param_provider import LigandBackedParameterisationProvider
from .mmff_provider import (
    _cx_residue_to_rdmol,
    _KCAL_TO_KJ, _ANG_TO_NM, _DEG_TO_RAD, _MDYNE_A_TO_KCAL,
    _ELEMENT_LJ, _DEFAULT_LJ,
)


_dgl_guarded = False


def _ensure_dgl_importable():
    '''
    Make ``import dgl`` succeed even when no ``graphbolt`` native library
    matches the installed PyTorch.

    ``graphbolt`` is DGL's optional GPU sampling subsystem; espaloma-charge
    uses only core ``DGLGraph``, never graphbolt.  But DGL imports it eagerly
    (``dgl`` -> ``dataloading`` -> ``distributed`` -> ``dist_graph`` ->
    ``from .. import graphbolt as gb``), so a missing/incompatible graphbolt
    build takes the whole package down.

    If a plain ``import dgl`` fails, we register an empty stub module under
    ``sys.modules['dgl.graphbolt']``.  CPython's ``IMPORT_FROM`` falls back to
    a ``sys.modules`` lookup when the attribute is absent, so DGL's internal
    ``from .. import graphbolt`` resolves to the stub and the real (broken)
    ``graphbolt/__init__.py`` never executes.  Every ``gb.<name>`` reference in
    DGL's ``distributed`` package lives inside a function body, so nothing at
    import time needs the stub to be populated.

    Idempotent, and a no-op when DGL imports cleanly (patched install or a
    future torch-compatible DGL) so working setups are never disturbed.
    '''
    global _dgl_guarded
    if _dgl_guarded:
        return
    _dgl_guarded = True

    import sys
    if 'dgl' in sys.modules:
        return  # already imported successfully by someone else
    try:
        import dgl  # noqa: F401
        return     # clean import — nothing to do
    except Exception:
        # A failed import leaves no partial 'dgl' in sys.modules (CPython
        # removes it on error), so we can safely stub and let the caller's
        # own `import espaloma_charge`/`import dgl` retry against the stub.
        import types
        sys.modules.setdefault('dgl.graphbolt', types.ModuleType('dgl.graphbolt'))


def _esc_system_to_ffxml(residue, mol, charges) -> str:
    '''
    Build a ForceField XML template using the supplied *charges*, MMFF94
    bonds/angles/torsions, and the GAFF2-element LJ table.

    Atom type prefix is ``ESC_`` to prevent clashes with the ``MMF_`` types
    that ``amber14+mmff`` registers separately.
    '''
    import xml.etree.ElementTree as ET
    from rdkit.Chem import rdForceFieldHelpers

    cx_atoms  = list(residue.atoms)
    n_atoms   = mol.GetNumAtoms()
    rname     = residue.name
    type_names = [f'ESC_{rname}_{i}' for i in range(n_atoms)]

    props = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94')
    if props is None:
        raise RuntimeError(
            f'MMFF94 could not assign bond/angle parameters for residue {rname!r}')

    def _lj(i):
        return _ELEMENT_LJ.get(mol.GetAtomWithIdx(i).GetAtomicNum(), _DEFAULT_LJ)

    sigma   = [_lj(i)[0] for i in range(n_atoms)]
    epsilon = [_lj(i)[1] for i in range(n_atoms)]

    root = ET.Element('ForceField')

    at_el = ET.SubElement(root, 'AtomTypes')
    for i, cx_a in enumerate(cx_atoms):
        ET.SubElement(at_el, 'Type', {
            'name':    type_names[i],
            'class':   type_names[i],
            'element': cx_a.element.name.capitalize(),
            'mass':    str(cx_a.element.mass),
        })

    res_el = ET.SubElement(ET.SubElement(root, 'Residues'), 'Residue', name=rname)
    for i, cx_a in enumerate(cx_atoms):
        ET.SubElement(res_el, 'Atom',
            name=cx_a.name, type=type_names[i], charge=str(charges[i]))
    bonds_obj = residue.atoms.intra_bonds
    a1_col, a2_col = bonds_obj.atoms
    for a1_name, a2_name in zip(a1_col.names, a2_col.names):
        ET.SubElement(res_el, 'Bond', atomName1=a1_name, atomName2=a2_name)

    nb_el = ET.SubElement(root, 'NonbondedForce',
        coulomb14scale='0.8333333333333333', lj14scale='0.5')
    ET.SubElement(nb_el, 'UseAttributeFromResidue', name='charge')
    for i in range(n_atoms):
        ET.SubElement(nb_el, 'Atom',
            type=type_names[i], sigma=str(sigma[i]), epsilon=str(epsilon[i]))

    # HarmonicBondForce — MMFF94 bond stretches
    # GetMMFFBondStretchParams -> (bondType, kb, r0); kb in md/Å, r0 in Å.
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

    # HarmonicAngleForce — MMFF94 angle bends
    # GetMMFFAngleBendParams -> (angleType, ka, theta0); ka in md·Å/rad², deg.
    af_el = None
    for j in range(n_atoms):
        nbrs = [b.GetOtherAtomIdx(j) for b in mol.GetAtomWithIdx(j).GetBonds()]
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

    # PeriodicTorsionForce — MMFF94 proper torsions
    # MMFF94: E = (V1/2)*(1+cos φ) + (V2/2)*(1−cos 2φ) + (V3/2)*(1+cos 3φ)
    # Mapped to OpenMM k*(1+cos(nφ−phase)):
    #   V1 → k=|V1|/2, n=1, phase=0   (sign flip → phase += π)
    #   V2 → k=|V2|/2, n=2, phase=π
    #   V3 → k=|V3|/2, n=3, phase=0
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
                tkey = tuple(sorted([i, j, k, l]) + [j, k])
                if tkey in seen_torsions:
                    continue
                seen_torsions.add(tkey)
                tp = props.GetMMFFTorsionParams(mol, i, j, k, l)
                if tp is None:
                    continue
                # GetMMFFTorsionParams -> (torsionType, V1, V2, V3); V in kcal/mol.
                v1, v2, v3 = tp[1], tp[2], tp[3]
                if tf_el is None:
                    tf_el = ET.SubElement(root, 'PeriodicTorsionForce')
                for v, n_per, base_phase in ((v1, 1, 0.0), (v2, 2, math.pi), (v3, 3, 0.0)):
                    if abs(v) < 1e-6:
                        continue
                    k_kJ  = abs(v) * _KCAL_TO_KJ / 2.0
                    phase = base_phase + (math.pi if v < 0.0 else 0.0)
                    ET.SubElement(tf_el, 'Proper',
                        type1=type_names[i], type2=type_names[j],
                        type3=type_names[k], type4=type_names[l],
                        periodicity1=str(n_per),
                        phase1=str(phase),
                        k1=str(k_kJ))

    return ET.tostring(root, encoding='unicode')


# In-memory cache of the loaded GNN, shared across all residues in a session.
_ESPALOMA_MODEL = None


def _load_espaloma_model():
    '''
    Load (and cache) the espaloma-charge GNN.

    This reimplements the model-loading half of ``espaloma_charge.app.charge``
    to fix two shortcomings of the stock function for in-session use:

    - The stock ``charge()`` downloads/reads ``model.pt`` from a *relative*
      path in the current working directory, which is unpredictable inside a
      ChimeraX session.  We cache it under ``torch.hub`` instead — a stable,
      always-writable per-user location.
    - The stock ``charge()`` reloads the model from disk on every call; since
      ISOLDE parameterises one residue at a time, that would re-deserialise the
      GNN for every ligand.  We load it once and keep it in memory.

    PyTorch >= 2.6 defaults ``torch.load`` to ``weights_only=True`` and refuses
    to unpickle the published ``model.pt`` (a full ``nn.Module``, not a
    state-dict).  The model is from choderalab's trusted release, so
    ``weights_only=False`` is appropriate here.
    '''
    global _ESPALOMA_MODEL
    if _ESPALOMA_MODEL is not None:
        return _ESPALOMA_MODEL

    import os
    import torch
    from urllib import request
    _ensure_dgl_importable()
    from espaloma_charge.app import MODEL_URL

    cache_dir = os.path.join(torch.hub.get_dir(), 'espaloma_charge')
    os.makedirs(cache_dir, exist_ok=True)
    model_path = os.path.join(cache_dir, 'model.pt')
    if not os.path.exists(model_path):
        request.urlretrieve(MODEL_URL.strip(), model_path)

    _ESPALOMA_MODEL = torch.load(model_path, weights_only=False)
    return _ESPALOMA_MODEL


def _espaloma_charges(mol):
    '''
    Run the espaloma-charge GNN on *mol* and return an ``(n_atoms,)`` numpy
    array of partial charges.  Mirrors ``espaloma_charge.app.charge`` but uses
    the cached model from :func:`_load_espaloma_model`.
    '''
    _ensure_dgl_importable()
    from espaloma_charge.utils import from_rdkit_mol

    model = _load_espaloma_model()
    graph = from_rdkit_mol(mol)
    graph = model(graph)
    return graph.ndata['q'].cpu().detach().flatten().numpy()


def _espaloma_charge_parameterise_residue(residue, logger) -> str:
    '''
    Build the RDKit mol, obtain GNN partial charges via espaloma-charge,
    then assemble the ForceField XML using MMFF94 bonded parameters.
    '''
    from rdkit.Chem import AllChem

    mol = _cx_residue_to_rdmol(residue)
    result = AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    if result == -1:
        raise RuntimeError(
            f'MMFF94 setup failed for {residue.name!r} (needed for bonded terms). '
            f'The molecule may contain elements or bonding patterns not supported '
            f'by MMFF94.')

    charges = _espaloma_charges(mol).tolist()
    return _esc_system_to_ffxml(residue, mol, charges)


class EspalomaChargeParameterisationProvider(LigandBackedParameterisationProvider):
    '''
    Parameterisation provider using AMBER14 for standard residues and a
    espaloma-charge / MMFF94 hybrid for ligands.

    - Partial charges: espaloma-charge GNN (AM1-BCC quality, pip-installable)
    - Bonded terms:    MMFF94 via RDKit (already bundled with ChimeraX)
    - LJ parameters:  GAFF2-element table, compatible with AMBER14 mixing

    Install the extra dependency with::

        pip install espaloma-charge dgl
    '''

    @property
    def backend_forcefield_key(self):
        return 'amber14+espaloma-charge'

    def _parameterise_one(self, residue, logger):
        xml = _espaloma_charge_parameterise_residue(residue, logger)
        return xml, {}
