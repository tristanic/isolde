# @Author: Tristan Croll
# @Date:   19-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 19-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Characterise ISOLDE's *current* forcefield-template resolution -- the path the
future ML parameterisation backend will replace. Driven through the read-only
``isolde preflight parameters`` dry-run, which resolves templates without ever
building an OpenMM Context.
'''

import pytest

pytestmark = pytest.mark.fast


def _preflight(session, model):
    from chimerax.isolde.validation.cmd import isolde_preflight_parameters
    # forcefield=None -> ISOLDE's configured default (defaults.OPENMM_FORCEFIELD)
    return isolde_preflight_parameters(session, model=model, forcefield=None)


def _residue_key(chain_id, name, number):
    return '/{} {} {}'.format(chain_id, name, number)


def test_preflight_schema_and_counts(session, model, golden):
    '''The dry-run dict's shape and per-class residue counts for the fixture.'''
    result = _preflight(session, model)

    for key in ('model', 'forcefield', 'n_total', 'n_matched', 'n_ambiguous',
                'n_unmatched', 'ready_for_simulation', 'unmatched', 'ambiguous'):
        assert key in result, 'preflight result lost key {!r}'.format(key)

    assert isinstance(result['ready_for_simulation'], bool)
    assert result['n_total'] == len(model.residues)
    assert (result['n_matched'] + result['n_ambiguous'] + result['n_unmatched']
            == result['n_total'])

    golden('template_resolution_counts', {
        'n_total': result['n_total'],
        'n_matched': result['n_matched'],
        'n_ambiguous': result['n_ambiguous'],
        'n_unmatched': result['n_unmatched'],
        'ready_for_simulation': result['ready_for_simulation'],
    })


def test_per_residue_template_map(session, model, golden):
    '''The explicit residue->template assignment is the heart of current
    parameterisation; pin it residue by residue.'''
    from chimerax.atomic import Residues
    from chimerax.isolde.validation.cmd import _get_forcefield
    from chimerax.isolde.openmm.openmm_interface import find_residue_templates

    ff, ligand_db, ff_name = _get_forcefield(session, None)
    residues = Residues(sorted(model.residues,
        key=lambda r: (r.chain_id, r.number, r.insertion_code)))
    tdict = find_residue_templates(residues, ff, ligand_db=ligand_db,
        logger=session.logger)

    mapping = {}
    for i, tname in tdict.items():
        r = residues[i]
        mapping[_residue_key(r.chain_id, r.name, r.number)] = tname
    golden('template_resolution_map', mapping)


def test_unmatched_residue_is_flagged(session, model):
    '''Today OpenMM's template system is the safety net: an unrecognisable
    residue is reported unmatched (the rewrite removes this net). Forcing a
    heavy atom to boron guarantees no template can match by name or topology.'''
    from chimerax.atomic import Element

    target = None
    for r in model.residues:
        heavy = r.atoms[r.atoms.element_names != 'H']
        if len(heavy):
            target = (r, heavy[0])
            break
    assert target is not None, 'fixture has no heavy atoms?'
    r, atom = target
    atom.element = Element.get_element('B')

    result = _preflight(session, model)
    assert result['ready_for_simulation'] is False

    hit = [u for u in result['unmatched']
           if u['number'] == r.number and u['chain_id'] == r.chain_id]
    assert hit, 'boron-mutated residue {} {}{} was not flagged unmatched'.format(
        r.name, r.chain_id, r.number)
    # The dry-run still offers actionable candidates keyed on the residue name.
    assert isinstance(hit[0]['candidates_by_name'], list)
