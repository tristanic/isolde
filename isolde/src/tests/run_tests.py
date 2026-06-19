# @Author: Tristan Croll
# @Date:   19-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 19-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Entry point for ISOLDE's characterisation test suite. Run inside ChimeraX:

    ChimeraX --nogui --exit --cmd "runscript /path/to/isolde/src/tests/run_tests.py [pytest args]"

Examples:
    ... run_tests.py -m fast          # only the no-OpenMM tests
    ... run_tests.py                  # the full suite
    ... run_tests.py --regenerate     # (re)write golden snapshots, then skip compares

ChimeraX's ``runscript`` injects ``session`` into this script's globals; we
stash it for the pytest fixtures, then hand the remaining argv to pytest.
'''

import os
import sys

here = os.path.dirname(os.path.abspath(__file__))

# Hand the live session to conftest's fixtures.
from chimerax.isolde.tests import _runtime
_runtime.session = session       # noqa: F821 -- provided by ChimeraX runscript

args = list(sys.argv[1:])
if '--regenerate' in args:
    args.remove('--regenerate')
    os.environ['ISOLDE_TEST_REGENERATE'] = '1'

# Always confine collection to this directory and pin the rootdir, so pytest
# can't wander up into the repo (e.g. extern/pybind11/tests) and so our
# pytest.ini (markers, import-mode) is the one that loads. Extra args from the
# caller (-m fast, -k ..., -v) are appended and still apply.
args = [here, '--rootdir', here] + args

try:
    import pytest
except ImportError:
    session.logger.error(                       # noqa: F821
        "pytest is not installed in ChimeraX's Python. Install it with:\n"
        "    ChimeraX -m pip install pytest")
    raise SystemExit(2)

exit_code = int(pytest.main(args))
session.logger.info('ISOLDE test suite finished with exit code {}'.format(exit_code))  # noqa: F821
# Propagate only failures, so a clean run doesn't surface a spurious SystemExit.
if exit_code != 0:
    raise SystemExit(exit_code)
