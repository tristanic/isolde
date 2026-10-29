# @Author: Tristan Croll
# @Date:   29-Oct-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 29-Oct-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Unit test for the symmetry-shell operator classification
(``symmetry_sim._operator_has_fixed_point``): only an operator instance that
leaves some point fixed can map a residue onto itself, so only those are tested
for special-position self-overlap.

Protein space groups contain only rotations, screws and translations, but
centrosymmetric small-molecule groups add improper operators (inversion, mirror,
glide, rotoinversion), which have no rotation axis. The old proper-rotation test
crashed on them (UnboundLocalError in Clipper's axis decomposition; COD 2226487,
P-1). Inversion centres and mirrors are exactly where small molecules sit on
special positions, so they must count as fixed-point operators.

Also builds a real symmetry shell for Clipper's bundled C2/c test crystal
(inversion, c-glide and a 2-fold, with a Cu on the 2-fold). Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_symmetry_fixed_point.py

Prints PASS/FAIL and exits non-zero on failure.
'''
import os
import numpy


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def _op(r, t):
    m = numpy.zeros((3, 4))
    m[:, :3] = r
    m[:, 3] = t
    return m


def _rot_z(deg):
    c, s = numpy.cos(numpy.radians(deg)), numpy.sin(numpy.radians(deg))
    return numpy.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def test_classification():
    from chimerax.isolde.openmm.symmetry_sim import _operator_has_fixed_point
    I = numpy.identity(3)
    two_z = numpy.diag([-1., -1., 1.])
    mirror_z = numpy.diag([1., 1., -1.])
    three_111 = numpy.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])
    cases = [
        ('pure translation', _op(I, (5, 0, 0)), False),
        ('2-fold through origin', _op(two_z, (0, 0, 0)), True),
        # x -> Rx + (I-R)p rotates about an axis through p.
        ('2-fold through (3,4,0)', _op(two_z, (I - two_z) @ (3, 4, 0)), True),
        ('2_1 screw', _op(two_z, (0, 0, 5)), False),
        ('2-fold + lattice shift along axis', _op(two_z, (6, 8, 10)), False),
        ('3-fold about [111]', _op(three_111, (0, 0, 0)), True),
        ('3_1 screw about [111]', _op(three_111, (2, 2, 2)), False),
        ('6-fold, off-origin axis', _op(_rot_z(60), (1, 0, 0)), True),
        ('inversion at origin', _op(-I, (0, 0, 0)), True),
        ('inversion at (0.5,1,1.5)', _op(-I, (1, 2, 3)), True),
        ('mirror at z=2', _op(mirror_z, (0, 0, 4)), True),
        ('glide in z=2', _op(mirror_z, (3, 0, 4)), False),
        ('-4 rotoinversion', _op(-_rot_z(90), (1, 1, 1)), True),
        ('-3 rotoinversion', _op(-three_111, (1, -2, 3)), True),
        ('-6 rotoinversion', _op(-_rot_z(60), (2, 0, 7)), True),
    ]
    for name, m, expected in cases:
        got = _operator_has_fixed_point(m)
        if got != expected:
            _fail('%s: _operator_has_fixed_point gave %s, expected %s' % (name, got, expected))
        print('  ok  %-34s -> %s' % (name, got))
    print('PASS test_classification')


def test_centrosymmetric_shell(session):
    from chimerax.isolde.openmm.symmetry_sim import (build_symmetry_shell,
        symmetry_sims_available)
    if not symmetry_sims_available():
        print('SKIP test_centrosymmetric_shell: OpenMM has no SymmetrySite')
        return
    import chimerax.clipper
    from chimerax.clipper.io.small_molecule import open_small_molecule_cif
    cif = os.path.join(os.path.dirname(chimerax.clipper.__file__), 'tests', 'cod_1100908.cif')
    model = open_small_molecule_cif(session, cif)
    try:
        shell = build_symmetry_shell(model, model.atoms, 8.0)
        if shell is None or not len(shell.copies):
            _fail('no symmetry copies sourced for C2/c test crystal')
        improper = {c.symop_index for c in shell.copies
                    if numpy.linalg.det(shell.symmats[c.symop_index][:, :3]) < 0}
        if not improper:
            _fail('no copies under an improper operator - test is not exercising anything')
        print('PASS test_centrosymmetric_shell: %d copies, %d under improper operators'
              % (len(shell.copies), sum(1 for c in shell.copies if c.symop_index in improper)))
    finally:
        session.models.close([model])


test_classification()
if 'session' in globals():
    test_centrosymmetric_shell(session)
