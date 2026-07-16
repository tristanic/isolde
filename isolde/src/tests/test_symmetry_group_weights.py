# @Author: Tristan Croll
# @Date:   16-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 16-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Regression test for the crystallographic symmetry group-weight table
(``symmetry_group_weight_table`` in ``openmm.symmetry_sim``), which weights the
nonbonded / GBSA copy<->copy interactions so each unique crystal contact is
counted exactly once. Pure arithmetic -- no GUI, no OpenMM context. Run inside
ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_symmetry_group_weights.py

Prints PASS/FAIL and exits non-zero on failure.
'''
import numpy


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def _trans(x, y, z):
    '''Pure-translation operator (R = I) as a (3,4) [R|t].'''
    m = numpy.zeros((3, 4))
    m[:3, :3] = numpy.eye(3)
    m[:, 3] = (x, y, z)
    return m


def _unflat(table, n):
    '''Recover w[i][j] from the flattened row-major table (index i + n*j).'''
    return [[table[i + n * j] for j in range(n)] for i in range(n)]


def run(session=None):
    from chimerax.isolde.openmm.symmetry_sim import (
        symmetry_group_weight_table, _op_inverse, _op_compose, _op_equal)

    def W(ops):
        return _unflat(symmetry_group_weight_table(ops), len(ops))

    def classes_sum_to_one(ops, tol=1e-6):
        '''Every relationship class {S,S^-1} present among all group-pairs (incl.
        the real legs) must have its weights sum to 1 -- the correctness
        invariant.'''
        n = len(ops)
        w = W(ops)
        pairs = [(x, y) for x in range(n) for y in range(x + 1, n)]
        S = {(x, y): _op_compose(_op_inverse(ops[x]), ops[y]) for (x, y) in pairs}
        seen = []
        for p in pairs:
            if any(_op_equal(S[p], s, 1e-3) or _op_equal(S[p], _op_inverse(s), 1e-3)
                   for s in seen):
                continue
            seen.append(S[p])
            members = [q for q in pairs
                       if _op_equal(S[q], S[p], 1e-3)
                       or _op_equal(S[q], _op_inverse(S[p]), 1e-3)]
            if abs(sum(w[x][y] for (x, y) in members) - 1.0) > tol:
                return False
        return True

    # (a) symmetric 2-body: a 2-fold operator (self-inverse) -- the validated case.
    twofold = numpy.zeros((3, 4)); twofold[:3, :3] = numpy.diag([-1.0, -1.0, 1.0])
    w = W([_trans(0, 0, 0), twofold])
    if not (abs(w[0][0] - 1.0) < 1e-9 and abs(w[1][1]) < 1e-9
            and abs(w[0][1] - 0.5) < 1e-9):
        _fail('2-fold {I,P}: expected real-real 1, real-copy 1/2, copy-copy 0; got %r' % w)
    print('PASS: 2-fold {I,P} -> real-real 1, real-copy 1/2, copy-copy(self) 0')

    # (b) whole-construct 3-way: linking operator S = P^-1 Q present -> copy-copy 0.
    I = _trans(0, 0, 0); P = _trans(10, 0, 0); Q = _trans(0, 10, 0)
    S = _op_compose(_op_inverse(P), Q)
    ops_b = [I, P, Q, S, _op_inverse(S), _op_inverse(P), _op_inverse(Q)]
    w = W(ops_b)
    if abs(w[1][2]) >= 1e-9:
        _fail('whole-construct 3-way: copy_P<->copy_Q should be 0, got %r' % w[1][2])
    if not classes_sum_to_one(ops_b):
        _fail('whole-construct 3-way: a relationship class does not sum to 1')
    print('PASS: whole-construct 3-way (P^-1Q present) -> copy_P<->copy_Q = 0, classes sum 1')

    # (c) local 3-way: P,Q,P^-1,Q^-1 present but NOT P^-1Q -- the bug case.
    ops_c = [I, P, Q, _op_inverse(P), _op_inverse(Q)]
    w = W(ops_c)
    if w[1][2] <= 1e-9:
        _fail('local 3-way: copy_P<->copy_Q must now be counted (>0), got %r' % w[1][2])
    if not classes_sum_to_one(ops_c):
        _fail('local 3-way: a relationship class does not sum to 1')
    print('PASS: local 3-way (P^-1Q absent) -> copy_P<->copy_Q = %.3f (>0), classes sum 1'
          % w[1][2])

    # (d) minimal local {I,P,Q} -> copy_P<->copy_Q is the sole rep -> weight 1.
    w = W([I, P, Q])
    if abs(w[1][2] - 1.0) >= 1e-9:
        _fail('minimal local {I,P,Q}: sole copy-copy rep must be weight 1, got %r' % w[1][2])
    print('PASS: minimal local {I,P,Q} -> copy_P<->copy_Q = 1 (sole representation)')

    # invariants across all sets: symmetric, in [0,1], diagonal 1 (real) / 0 (copy).
    for ops in (ops_b, ops_c, [I, P, Q], [_trans(0, 0, 0), twofold]):
        n = len(ops); w = W(ops)
        if abs(w[0][0] - 1.0) >= 1e-9:
            _fail('real-real weight != 1')
        for k in range(1, n):
            if abs(w[k][k]) >= 1e-9:
                _fail('copy self-weight != 0')
        for i in range(n):
            for j in range(n):
                if abs(w[i][j] - w[j][i]) >= 1e-12:
                    _fail('weight table not symmetric')
                if not (-1e-12 <= w[i][j] <= 1.0 + 1e-9):
                    _fail('weight %r out of [0,1]' % w[i][j])
    print('PASS: invariants (symmetric, [0,1], diagonal 1/0) hold')

    print('ALL PASS')


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
