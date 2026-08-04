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
  * the operator of its most-interior drawn image (nearest the ROI centroid).

Because every crystallographic representation of an atom folds the identical
force back to it, one full-coupling term is exact -- there is no fractional
``/n_sym`` splitting. And because a symmetry shell only exists for
crystallographic data (map defined throughout the cell) and is built
whole-residue, every mobile atom must have a covered representation -- a mobile
atom with none is a coverage/shell bug and must raise.

Drives the real method with lightweight fakes (real atoms for
coordinates/indexing; a recording fake force), so it needs no running simulation
or GUI. Run inside ChimeraX:

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


# Operators: 0 = identity (unused), 1 = small translation (image lands near the
# ROI), 2 = large translation (image lands far). Shared by the scenarios.
def _symmats():
    s = numpy.zeros((3, 3, 4))
    for k in range(3):
        s[k, :3, :3] = numpy.eye(3)
    s[1, :, 3] = (5.0, 0.0, 0.0)
    s[2, :, 3] = (50.0, 0.0, 0.0)
    return s


def _run_selection(construct, coverage, copies):
    '''Call the real _add_mdff_symmetry_terms with fakes; MDFF atoms = the whole
    construct at full coupling 3.0. Returns (force, mdff_atoms, fake_handler).'''
    from chimerax.isolde.openmm.openmm_interface import SimHandler
    shell = SimpleNamespace(symmats=_symmats(), copies=copies)
    sh = SimpleNamespace(
        _atoms=construct,
        _sim_construct=SimpleNamespace(mobile_atoms=construct,
                                       symmetry_coverage_atoms=coverage),
        _mdff_symmetry_term_indices={})
    idx = construct.indices(construct)
    ks = numpy.full(len(idx), 3.0)
    enableds = numpy.ones(len(idx))
    mdff_atoms = SimpleNamespace()
    f = _RecordingForce()
    SimHandler._add_mdff_symmetry_terms(sh, f, 'vol', shell, mdff_atoms,
                                        idx, ks, enableds)
    return f, mdff_atoms, sh


def run(session):
    from chimerax.core.commands import run as run_cmd
    from chimerax.isolde.openmm.openmm_interface import _symmat_to_transform12
    from chimerax.isolde.openmm.custom_forces import SymmetryAwareCubicInterpMapForce

    m = run_cmd(session, 'open "%s"' % FIXTURE.replace('\\', '/'), log=False)[0]
    try:
        atoms = m.atoms
        Copy = lambda atom, op: SimpleNamespace(parent_atom=atom, symop_index=op)

        # === (A) normal: 12 covered ROI atoms + 3 distant parents, all with a
        #         covered representation. p14 has TWO images (op 2 listed FIRST)
        #         so we distinguish "pick nearest ROI" from "pick first". ===
        construct = atoms[:15]
        coverage = construct[:12]
        p12, p13, p14 = construct[12], construct[13], construct[14]
        copies = [Copy(p12, 1), Copy(p14, 2), Copy(p14, 1), Copy(p13, 1)]
        f, mdff_atoms, sh = _run_selection(construct, coverage, copies)

        if len(f.calls) != 1:
            _fail('expected exactly one add_atoms call, got %d' % len(f.calls))
        call = f.calls[0]
        termed = list(call['indices'])

        if sorted(termed) != list(range(15)):
            _fail('termed atoms %r != 0..14' % sorted(termed))
        if len(termed) != len(set(termed)):
            _fail('an atom received more than one term: %r' % termed)
        print('PASS: exactly one term per atom (covered ROI + represented parent)')

        if not numpy.allclose(call['ks'], 3.0):
            _fail('couplings are not full (3.0): %r' % call['ks'])
        print('PASS: every term carries the full coupling constant (no /n_sym split)')

        tf_by_atom = {int(i): call['transforms'][row]
                      for row, i in enumerate(termed)}
        ident = SymmetryAwareCubicInterpMapForce.IDENTITY_TRANSFORM
        for i in range(12):
            if not numpy.allclose(tf_by_atom[i], ident):
                _fail('covered atom %d did not get the identity transform' % i)
        near = _symmat_to_transform12(_symmats()[1])
        for i in (12, 13, 14):
            if not numpy.allclose(tf_by_atom[i], near):
                _fail('parent %d did not get the nearest-image (op 1) transform' % i)
        print('PASS: covered atoms use identity; distant parents use their '
              'most-interior image operator (op 1, not the closer-listed op 2)')

        si = numpy.asarray(mdff_atoms.sim_indices)
        if (si < 0).any():
            _fail('some represented atom has sim_index -1: %r' % list(si))
        d = sh._mdff_symmetry_term_indices['vol']
        if set(d) != set(termed):
            _fail('term map keys %r != termed atoms %r' % (sorted(d), sorted(termed)))
        for i, entry in d.items():
            if not (isinstance(entry, tuple) and len(entry) == 2):
                _fail('term map entry for %d is not a single (fi, tf): %r' % (i, entry))
        print('PASS: sim_indices + term map record exactly one term per atom')

        # === (B) invariant guard: a mobile atom that is neither covered nor has
        #         any image must raise (should never happen in a correct
        #         crystallographic sim). Atoms 12-14 keep their images; atom 15
        #         is non-covered with no image. ===
        construct = atoms[:16]
        coverage = construct[:12]
        copies = [Copy(construct[12], 1), Copy(construct[13], 1),
                  Copy(construct[14], 1)]
        raised = False
        try:
            _run_selection(construct, coverage, copies)
        except RuntimeError:
            raised = True
        if not raised:
            _fail('a mobile atom with no covered representation did not raise')
        print('PASS: a mobile atom with no covered representation raises '
              '(crystallographic invariant guard)')

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
