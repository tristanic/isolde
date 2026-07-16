# @Author: Tristan Croll
# @Date:   16-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 16-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Unit test for the symmetry-aware MDFF term-SELECTION logic
(``SimHandler._add_mdff_symmetry_terms``): each real atom must get EXACTLY ONE
map term at its full coupling constant, with the transform chosen as

  * identity, when the atom's own position is inside the covered map box, else
  * the operator of its most-interior drawn image (nearest the ROI centroid),
  * and no term at all when it has no covered representation.

Because every crystallographic representation of an atom folds the identical
force back to it, one full-coupling term is exact -- there is no fractional
``/n_sym`` splitting. This drives the real method with lightweight fakes (real
atoms for coordinates/indexing; a recording fake force), so it needs no running
simulation or GUI. Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_symmetry_mdff_selection.py

Prints PASS/FAIL and exits non-zero on failure.
'''
import os
import numpy
from types import SimpleNamespace

FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '1pmx_1.pdb')


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


class _RecordingForce:
    '''Stands in for a SymmetryAwareCubicInterpMapForce; records the one
    add_atoms call and hands back sequential force indices.'''
    def __init__(self):
        self.calls = []

    def add_atoms(self, indices, ks, enableds, transforms=None):
        self.calls.append(dict(indices=numpy.asarray(indices),
                               ks=numpy.asarray(ks),
                               enableds=numpy.asarray(enableds),
                               transforms=numpy.asarray(transforms)))
        return numpy.arange(len(indices), dtype=numpy.int32)


def run(session):
    from chimerax.core.commands import run as run_cmd
    from chimerax.isolde.openmm.openmm_interface import (
        SimHandler, _symmat_to_transform12)
    from chimerax.isolde.openmm.custom_forces import SymmetryAwareCubicInterpMapForce

    m = run_cmd(session, 'open "%s"' % FIXTURE.replace('\\', '/'), log=False)[0]
    try:
        all_atoms = m.atoms[:20]              # the "construct"
        coverage = all_atoms[:12]             # region of interest (own pos covered)
        # Atoms 12,13,14 are non-covered "distant parents" with drawn images;
        # atoms 15-19 are non-covered with NO image (should get no term).
        p12, p13, p14 = all_atoms[12], all_atoms[13], all_atoms[14]

        # Operators: 0 = identity (unused here), 1 = small translation (image lands
        # near the ROI), 2 = large translation (image lands far).
        symmats = numpy.zeros((3, 3, 4))
        for k in range(3):
            symmats[k, :3, :3] = numpy.eye(3)
        symmats[1, :, 3] = (5.0, 0.0, 0.0)
        symmats[2, :, 3] = (50.0, 0.0, 0.0)

        # p14 has TWO images (op 2 listed FIRST) so the test distinguishes
        # "pick nearest ROI" from "pick first".
        Copy = lambda atom, op: SimpleNamespace(parent_atom=atom, symop_index=op)
        shell = SimpleNamespace(
            symmats=symmats,
            copies=[Copy(p12, 1), Copy(p14, 2), Copy(p14, 1), Copy(p13, 1)])

        sh = SimpleNamespace(
            _atoms=all_atoms,
            _sim_construct=SimpleNamespace(mobile_atoms=all_atoms,
                                           symmetry_coverage_atoms=coverage),
            _mdff_symmetry_term_indices={})

        indices = all_atoms.indices(all_atoms)                 # 0..19
        ks = numpy.full(len(all_atoms), 3.0)                   # full coupling
        enableds = numpy.ones(len(all_atoms))
        mdff_atoms = SimpleNamespace()
        f = _RecordingForce()

        SimHandler._add_mdff_symmetry_terms(
            sh, f, 'vol', shell, mdff_atoms, indices, ks, enableds)

        if len(f.calls) != 1:
            _fail('expected exactly one add_atoms call, got %d' % len(f.calls))
        call = f.calls[0]
        termed = list(call['indices'])

        # --- one term per covered ROI atom + per represented parent, none else ---
        expected = list(range(12)) + [12, 13, 14]     # 15..19 have no image
        if sorted(termed) != sorted(expected):
            _fail('termed atoms %r != expected %r' % (sorted(termed), expected))
        if len(termed) != len(set(termed)):
            _fail('an atom received more than one term: %r' % termed)
        print('PASS: exactly one term per covered atom / represented parent; '
              'atoms with no covered representation get none')

        # --- full coupling, never fractional ---
        if not numpy.allclose(call['ks'], 3.0):
            _fail('couplings are not full (3.0): %r' % call['ks'])
        print('PASS: every term carries the full coupling constant (no /n_sym split)')

        # --- transform choice ---
        tf_by_atom = {int(i): call['transforms'][row]
                      for row, i in enumerate(termed)}
        ident = SymmetryAwareCubicInterpMapForce.IDENTITY_TRANSFORM
        for i in range(12):
            if not numpy.allclose(tf_by_atom[i], ident):
                _fail('covered atom %d did not get the identity transform' % i)
        near = _symmat_to_transform12(symmats[1])
        for i in (12, 13, 14):
            if not numpy.allclose(tf_by_atom[i], near):
                _fail('parent %d did not get the nearest-image (op 1) transform' % i)
        print('PASS: covered atoms use identity; distant parents use their '
              'most-interior image operator (op 1, not the closer-listed op 2)')

        # --- sim_indices: termed atoms get a real index, the rest -1 ---
        si = numpy.asarray(mdff_atoms.sim_indices)
        for i in range(20):
            if i in termed and si[i] < 0:
                _fail('termed atom %d has sim_index -1' % i)
            if i not in termed and si[i] != -1:
                _fail('unrepresented atom %d has sim_index %d (should be -1)'
                      % (i, si[i]))
        # --- term map: one (force_index, transform) per termed atom ---
        d = sh._mdff_symmetry_term_indices['vol']
        if set(d) != set(termed):
            _fail('term map keys %r != termed atoms %r' % (sorted(d), sorted(termed)))
        for i, entry in d.items():
            if not (isinstance(entry, tuple) and len(entry) == 2):
                _fail('term map entry for %d is not a single (fi, tf): %r' % (i, entry))
        print('PASS: sim_indices + term map record exactly one term per atom')

        print('ALL PASS')
    finally:
        session.models.close([m])


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
