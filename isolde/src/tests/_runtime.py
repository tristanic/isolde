# @Author: Tristan Croll
# @Date:   19-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 19-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Tiny holder used to hand the live ChimeraX ``session`` from ``run_tests.py``
(executed via ``ChimeraX ... runscript``) to the pytest fixtures in
``conftest.py``. ChimeraX's bundled Python has no notion of a global session,
so we stash it here at startup and read it back from the ``session`` fixture.
'''

session = None
