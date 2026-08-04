# Headless pure-logic tests for settle_poses (no simulation / GUI / map required).
#
# Exercises the parts of settle_poses that are pure Python -- pose validation, the
# do-no-harm winner selection, and SettleResult assembly. The force-math invariant that
# E_pack recovers a pair interaction is covered by test_nb_group_softcore.py; the full
# end-to-end settle needs a GUI simulation and is covered by the GUI harness.

import numpy as np

from chimerax.isolde.refine.settle_poses import (
    _validate_poses, _choose_winner, SettleResult, CURRENT_LABEL)
from chimerax.core.errors import UserError

_failures = []


def _check(cond, msg):
    if cond:
        print("  ok:", msg)
    else:
        print("  FAIL:", msg)
        _failures.append(msg)


def _raises(fn, exc=UserError):
    try:
        fn()
    except exc:
        return True
    except Exception:
        return False
    return False


def _cand(e, label, is_current):
    # (E_pack, label, settled_construct_coords, is_current) -- coords unused by _choose_winner
    return (e, label, np.zeros((1, 3)), is_current)


def run(session=None):
    _dlog = lambda *a, **k: None

    # ---- _validate_poses ----------------------------------------------------
    n = 3
    good = _validate_poses([('a', np.zeros((n, 3))),
                            ('b', [[0, 0, 0], [1, 1, 1], [2, 2, 2]])], n, UserError)
    _check(len(good) == 2 and good[0][0] == 'a' and good[0][1].shape == (n, 3),
           "_validate_poses accepts (label, (n,3)) and coerces to float arrays")
    _check(_raises(lambda: _validate_poses([('bad', np.zeros((2, 3)))], n, UserError)),
           "_validate_poses rejects a wrong-shaped pose")
    _check(_raises(lambda: _validate_poses([42], n, UserError)),
           "_validate_poses rejects a non-(label, coords) entry")
    _check(_raises(lambda: _validate_poses([], n, UserError)),
           "_validate_poses rejects an empty pose list")

    # ---- _choose_winner: do-no-harm -----------------------------------------
    # best pose IS current -> keep current
    cands = [_cand(-5.0, CURRENT_LABEL, True), _cand(-3.0, 'A', False)]
    best, kept = _choose_winner(cands, 1.0, _dlog)
    _check(kept and best[1] == CURRENT_LABEL, "best==current -> keep current")

    # challenger beats current by >= margin -> commit challenger
    cands = [_cand(-10.0, 'A', False), _cand(-3.0, CURRENT_LABEL, True)]
    best, kept = _choose_winner(cands, 1.0, _dlog)
    _check((not kept) and best[1] == 'A',
           "challenger beats current by >= margin (7 >= 1) -> commit challenger")

    # challenger beats current by < margin -> keep current (do no harm)
    cands = [_cand(-3.5, 'A', False), _cand(-3.0, CURRENT_LABEL, True)]
    best, kept = _choose_winner(cands, 1.0, _dlog)
    _check(kept and best[1] == CURRENT_LABEL,
           "challenger beats current by < margin (0.5 < 1) -> keep current")

    # no current present (defensive): still commits the best
    cands = [_cand(-4.0, 'A', False), _cand(-2.0, 'B', False)]
    best, kept = _choose_winner(cands, 1.0, _dlog)
    _check((not kept) and best[1] == 'A', "no current candidate -> commit best")

    _check(_choose_winner([], 1.0, _dlog) == (None, False),
           "_choose_winner([]) -> (None, False)")

    # ---- SettleResult defaults ----------------------------------------------
    r = SettleResult(committed='A', kept_current=False, applied=True,
                     energies=[('A', -4.0), ('B', -2.0)], n_culled=1, accept_margin_kJ=2.0)
    _check(r.map_decoupled is True and r.score_mode == 'E_pack'
           and r.committed_coords is None and r.metrics is None,
           "SettleResult defaults (map_decoupled=True, score_mode='E_pack', ...)")
    _check(r.energies[0] == ('A', -4.0) and r.n_culled == 1,
           "SettleResult carries energies + n_culled")

    if _failures:
        print(f"\n{len(_failures)} CHECK(S) FAILED")
        raise SystemExit(1)
    print("\nALL PASS")


if 'session' in dict(globals()) and session is not None:
    run(session)
