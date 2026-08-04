# @Author: Tristan Croll
# @Date:   04-Aug-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 04-Aug-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Headless test: a coordinated metal must never be turned into a chiral-centre
definition. RDKit CIP-flags chlorophyll's central Mg (the chlorin's reduced ring D
breaks the four-N coordination symmetry), but that is a coordination geometry, not
an organic sp3 stereocentre -- a chirality restraint on it fights the metal-site
coordination and mangles the ring. ``chiral_definitions_from_ccd`` /
``reference_cip_codes`` must skip it (the metal analogue of the phosphate spurious-
centre skip) while keeping CLA's genuine ring/backbone carbon stereocentres.

    run_chimerax.bat --nogui --exit --script src/tests/test_chirality_metal.py

Needs the ChemComp CCD store populated (CLA); SKIPs cleanly if it is not.
'''


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def run(session):
    try:
        from rdkit import Chem
    except ImportError:
        print('SKIP: RDKit not available')
        return
    from chimerax.isolde.atomic import chirality as ch, rdkit_bridge as rb

    mol, status = rb.template_to_rdkit(session, 'CLA')
    if mol is None:
        print('SKIP: CLA CCD template unavailable (empty ChemComp store)')
        return
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    # Precondition: RDKit really does CIP-flag the Mg (else there is nothing to skip).
    mg_atoms = [a for a in mol.GetAtoms() if a.GetSymbol() == 'Mg']
    if not mg_atoms:
        _fail('CLA RDKit build has no Mg')
    if not any(a.HasProp('_CIPCode') for a in mg_atoms):
        print('SKIP: RDKit did not CIP-flag CLA Mg on this build; nothing to exclude')
        return

    # The generator must exclude the metal centre but keep the real carbon centres.
    defs = ch.chiral_definitions_from_ccd(session, 'CLA')
    if 'MG' in defs:
        _fail('CLA Mg was generated as a chiral centre (must be excluded as a metal)')
    non_metal = [d for d in defs if d != 'MG']
    if not non_metal:
        _fail('CLA lost all chiral centres (the real ring stereocentres must survive)')

    # The informational reference CIP codes must likewise skip the metal.
    codes = ch.reference_cip_codes(session, 'CLA')
    if 'MG' in codes:
        _fail('CLA Mg appears in reference CIP codes (must be excluded)')

    # Helper sanity: a metal is a metal centre, a carbon is not.
    c_atom = next(a for a in mol.GetAtoms() if a.GetSymbol() == 'C')
    if not ch._is_metal_centre(mg_atoms[0]):
        _fail('_is_metal_centre returned False for Mg')
    if ch._is_metal_centre(c_atom):
        _fail('_is_metal_centre returned True for carbon')

    print('PASS: coordinated metal (CLA Mg) excluded from chiral defs (%d ring centres kept)'
          % len(non_metal))
    print('ALL PASS')


try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
