# @Author: Tristan Croll
# @Date:   19-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 19-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Characterise the stable edges of ISOLDE's command surface: the read-only
``isolde status`` report and the custom ``IsoldeStructureArg`` argument parser
(single-structure selection + its error behaviour). These are the contracts
scripts and the GUI depend on, independent of the simulation internals.
'''

import pytest

pytestmark = pytest.mark.fast


def test_isolde_status_dict(session):
    '''Read-only status works even before "isolde start" and returns the
    documented shape.'''
    from chimerax.isolde.cmd.cmd import isolde_status
    result = isolde_status(session)
    for key in ('isolde_started', 'selected_model', 'simulation_running',
                'forcefield'):
        assert key in result
    assert isinstance(result['isolde_started'], bool)
    assert isinstance(result['simulation_running'], bool)


def test_isolde_structure_arg_parses_single_model(session, model):
    from chimerax.isolde.cmd.argspec import IsoldeStructureArg
    parsed, text, rest = IsoldeStructureArg.parse(
        '#{}'.format(model.id_string), session)
    assert parsed is model


def test_isolde_structure_arg_rejects_no_match(session, model):
    from chimerax.isolde.cmd.argspec import IsoldeStructureArg
    from chimerax.core.commands import AnnotationError
    # A spec that resolves to no atomic structure must be rejected, not guessed.
    with pytest.raises(AnnotationError):
        IsoldeStructureArg.parse('#999', session)
