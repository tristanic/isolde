# @Author: Claude
# @Date:   03-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 03-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Regression test for the MMFF94 parameter extraction shared by the
``amber14+mmff`` and ``amber14+espaloma-charge`` backends
(``mmff_provider.py`` / ``espaloma_charge_provider.py``).  Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_mmff_parameterisation.py

RDKit's ``MMFFMolProperties`` exposes the per-term getters as *methods on the
properties object* (``props.GetMMFFBondStretchParams(mol, i, j)``), returning the
*raw* MMFF94 parameters: ``(bondType, kb, r0)`` with ``kb`` in md/Å, and
``(angleType, ka, theta0)`` with ``ka`` in md·Å/rad².  Two failure modes this
guards against, both of which produce a valid-looking XML and a clean run while
being physically wrong:

  1. Calling the getters as free functions in ``rdForceFieldHelpers`` (they are
     not there) or reading the wrong tuple index (element 0 is the *type*, not
     the force constant).
  2. Omitting the md->kcal conversion constant (143.9325 kcal·mol⁻¹ per md·Å):
     every bond and angle then comes out ~144x too weak, so the ligand is
     effectively unbonded and drifts apart in simulation with no traceback.

The test builds ethanol as an in-memory ChimeraX structure (no fixture file),
runs it through the *real* provider code, and asserts (a) the C-C bond force
constant lands at AMBER scale (~256,000 kJ/mol/nm², not ~1,800), and (b) OpenMM
loads the generated ForceField XML and builds a system with finite energy.

Prints PASS/FAIL per check and exits non-zero on the first failure.  A pure
unit test (rdkit + openmm only): it does not need espaloma-charge/DGL installed,
because the espaloma-charge provider's bonded/units logic is exercised directly
via ``_esc_system_to_ffxml`` with dummy charges.
'''
import io
import math
import xml.etree.ElementTree as ET


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def _build_ethanol(session):
    '''Create an in-memory ethanol AtomicStructure (residue 'LIG') with explicit
    hydrogens and a sensible 3D conformation.  Atom order matches the RDKit mol
    that ``_cx_residue_to_rdmol`` will rebuild from it.'''
    from chimerax.atomic import AtomicStructure
    from chimerax.atomic.struct_edit import add_bond
    from rdkit import Chem
    from rdkit.Chem import AllChem

    rd = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    AllChem.EmbedMolecule(rd, randomSeed=1)
    AllChem.MMFFOptimizeMolecule(rd)
    conf = rd.GetConformer()

    s = AtomicStructure(session, name='ethanol-test')
    r = s.new_residue('LIG', 'A', 1)
    counts = {}
    atoms = []
    for i in range(rd.GetNumAtoms()):
        sym = rd.GetAtomWithIdx(i).GetSymbol()
        counts[sym] = counts.get(sym, 0) + 1
        a = s.new_atom('%s%d' % (sym, counts[sym]), sym)
        p = conf.GetAtomPosition(i)
        a.coord = (p.x, p.y, p.z)
        r.add_atom(a)
        atoms.append(a)
    for b in rd.GetBonds():
        add_bond(atoms[b.GetBeginAtomIdx()], atoms[b.GetEndAtomIdx()])
    return s, r


def _bond_k_between_elements(xml, e1, e2):
    '''Return the HarmonicBondForce k (kJ/mol/nm²) for the first bond whose two
    atom types are elements {e1, e2}, or None.'''
    root = ET.fromstring(xml)
    type_elem = {t.get('name'): t.get('element') for t in root.iter('Type')}
    bf = root.find('HarmonicBondForce')
    if bf is None:
        return None
    want = sorted([e1, e2])
    for bond in bf:
        pair = sorted([type_elem[bond.get('type1')], type_elem[bond.get('type2')]])
        if pair == want:
            return float(bond.get('k'))
    return None


def _openmm_energy(residue, xml):
    '''Load *xml* as an OpenMM ForceField, build a system for *residue*, and
    return its potential energy (kJ/mol) at the residue's current coordinates.'''
    from openmm import app, unit, Context, LangevinIntegrator, Vec3

    ff = app.ForceField(io.StringIO(xml))
    top = app.Topology()
    ch = top.addChain()
    ommres = top.addResidue(residue.name, ch)
    cx_atoms = list(residue.atoms)
    omm_atoms = []
    for a in cx_atoms:
        el = app.element.Element.getByAtomicNumber(a.element.number)
        omm_atoms.append(top.addAtom(a.name, el, ommres))
    idx = {a: i for i, a in enumerate(cx_atoms)}
    b0, b1 = residue.atoms.intra_bonds.atoms
    for x, y in zip(b0, b1):
        top.addBond(omm_atoms[idx[x]], omm_atoms[idx[y]])

    system = ff.createSystem(top, nonbondedMethod=app.NoCutoff, constraints=None)
    positions = [Vec3(float(a.coord[0]) * 0.1, float(a.coord[1]) * 0.1,
                      float(a.coord[2]) * 0.1) for a in cx_atoms]
    integ = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                               0.001 * unit.picosecond)
    ctx = Context(system, integ)
    ctx.setPositions(positions * unit.nanometer)
    return ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)


def _check_provider(session, label, make_xml):
    '''Shared assertions for one provider's generated XML.'''
    s, r = _build_ethanol(session)
    try:
        xml = make_xml(r)
        cc_k = _bond_k_between_elements(xml, 'C', 'C')
        if cc_k is None:
            _fail('%s: no C-C bond found in generated XML' % label)
        # Buggy (missing md->kcal factor) gives ~1,782; correct is ~256,000.
        if not (150_000 < cc_k < 400_000):
            _fail('%s: C-C bond k=%.0f kJ/mol/nm^2 outside physical range '
                  '[150k, 400k]; likely a units regression (missing the '
                  '143.9325 md->kcal factor gives ~1,782)' % (label, cc_k))
        energy = _openmm_energy(r, xml)
        if not math.isfinite(energy):
            _fail('%s: OpenMM system energy is not finite' % label)
        print('PASS: %s -> C-C bond k=%.0f kJ/mol/nm^2 (AMBER scale), '
              'OpenMM system built, energy=%.2f kJ/mol' % (label, cc_k, energy))
    finally:
        session.models.close([s])


def _mmff_xml(residue):
    from chimerax.isolde.openmm.mmff_provider import _mmff_parameterise_residue
    return _mmff_parameterise_residue(residue, None)


def _esc_xml(residue):
    # Exercise the espaloma-charge provider's bonded/units logic without needing
    # espaloma-charge/DGL: build the mol as the provider does, then call
    # _esc_system_to_ffxml directly with dummy (zero) charges.
    from rdkit.Chem import AllChem
    from chimerax.isolde.openmm.mmff_provider import _cx_residue_to_rdmol
    from chimerax.isolde.openmm.espaloma_charge_provider import _esc_system_to_ffxml
    mol = _cx_residue_to_rdmol(residue)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    charges = [0.0] * mol.GetNumAtoms()
    return _esc_system_to_ffxml(residue, mol, charges)


def run(session):
    _check_provider(session, 'amber14+mmff', _mmff_xml)
    _check_provider(session, 'amber14+espaloma-charge (bonded/units)', _esc_xml)
    print('ALL PASS')


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
