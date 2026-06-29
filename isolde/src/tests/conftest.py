# @Author: Tristan Croll
# @Date:   19-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 19-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Shared fixtures for ISOLDE's characterisation test suite.

These tests pin down what ISOLDE *does today*, so a future rewrite of
simulation startup (ML-based parameterisation + model-resident force-field
store) makes any silent behaviour change fail loudly. Everything runs inside
ChimeraX's bundled Python via ``run_tests.py``.
'''

import os
import json
import math

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
GOLDEN_DIR = os.path.join(HERE, 'golden')
FIXTURE_PDB = os.path.join(HERE, '1pmx_1.pdb')


def cx_path(path):
    '''ChimeraX command parser is happier with forward slashes on Windows.'''
    return path.replace('\\', '/')


def _regenerating():
    return os.environ.get('ISOLDE_TEST_REGENERATE', '') not in ('', '0', 'false', 'False')


# ---------------------------------------------------------------------------
# Core fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope='session')
def session():
    '''The live ChimeraX session, injected by run_tests.py at startup.'''
    import _runtime
    if _runtime.session is None:
        raise RuntimeError(
            'No ChimeraX session available. Run this suite via:\n'
            '  ChimeraX --nogui --exit --cmd "runscript .../tests/run_tests.py"'
        )
    return _runtime.session


@pytest.fixture
def model(session):
    '''Open the bundled 1pmx fixture (small, already protonated) and close it
    on teardown. A fresh copy per test so in-test mutations don't leak.'''
    from chimerax.core.commands import run
    m = run(session, 'open "{}"'.format(cx_path(FIXTURE_PDB)))[0]
    yield m
    try:
        session.models.close([m])
    except Exception:
        pass


@pytest.fixture(scope='session')
def cpu_platform():
    '''A deterministic, GPU-free OpenMM platform for single-point energies.'''
    import openmm
    for name in ('Reference', 'CPU'):
        try:
            return openmm.Platform.getPlatformByName(name)
        except Exception:
            continue
    raise RuntimeError('Neither Reference nor CPU OpenMM platform is available.')


# ---------------------------------------------------------------------------
# Golden-snapshot helper
# ---------------------------------------------------------------------------

def _json_default(obj):
    # numpy scalars / arrays -> native python
    try:
        import numpy
        if isinstance(obj, numpy.generic):
            return obj.item()
        if isinstance(obj, numpy.ndarray):
            return obj.tolist()
    except ImportError:
        pass
    raise TypeError('Cannot JSON-encode {!r}'.format(type(obj)))


def _normalise(data):
    '''Round-trip through JSON so tuples become lists etc., matching what was
    stored on disk; comparison is then like-for-like.'''
    return json.loads(json.dumps(data, default=_json_default, sort_keys=True))


def _assert_close(expected, actual, rtol, atol, where):
    if isinstance(expected, dict):
        assert isinstance(actual, dict), '{}: type changed'.format(where)
        assert set(expected) == set(actual), \
            '{}: keys changed (added {}, removed {})'.format(
                where, set(actual) - set(expected), set(expected) - set(actual))
        for k in expected:
            _assert_close(expected[k], actual[k], rtol, atol, '{}.{}'.format(where, k))
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(expected) == len(actual), \
            '{}: list length changed ({} -> {})'.format(where, len(expected), len(actual))
        for i, (e, a) in enumerate(zip(expected, actual)):
            _assert_close(e, a, rtol, atol, '{}[{}]'.format(where, i))
    elif isinstance(expected, bool) or expected is None:
        assert expected == actual, '{}: {!r} -> {!r}'.format(where, expected, actual)
    elif isinstance(expected, (int, float)):
        assert isinstance(actual, (int, float)), '{}: {!r} -> {!r}'.format(where, expected, actual)
        tol = atol + rtol * abs(expected)
        assert math.isclose(expected, actual, rel_tol=rtol, abs_tol=tol), \
            '{}: {!r} -> {!r} (tol {})'.format(where, expected, actual, tol)
    else:
        assert expected == actual, '{}: {!r} -> {!r}'.format(where, expected, actual)


@pytest.fixture
def golden():
    '''Returns ``check(name, data, rtol=, atol=)``. On the first run (or with
    --regenerate) it writes the snapshot and skips; thereafter it asserts the
    current output matches the committed snapshot within tolerance.'''
    def check(name, data, rtol=1e-3, atol=1e-6):
        data = _normalise(data)
        path = os.path.join(GOLDEN_DIR, name + '.json')
        if _regenerating() or not os.path.exists(path):
            os.makedirs(GOLDEN_DIR, exist_ok=True)
            with open(path, 'w') as f:
                json.dump(data, f, indent=2, sort_keys=True)
            pytest.skip('golden "{}" (re)generated -- rerun to compare'.format(name))
        with open(path) as f:
            expected = json.load(f)
        _assert_close(expected, data, rtol, atol, name)
    return check
