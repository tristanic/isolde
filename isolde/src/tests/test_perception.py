# @Author: Tristan Croll
# @Date:   19-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 19-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Characterise ISOLDE's existing chemical-perception helpers -- comparing modelled
residues against templates, and building a connectivity graph from a residue.
The future ML parameterisation layer will lean heavily on perception of exactly
this kind, so pin what's there now.
'''

import numpy
import pytest

pytestmark = pytest.mark.fast


def test_find_incorrect_residues(session, model, golden):
    '''Which residues today's template comparison flags as not matching their
    CCD/built-in template (heavy atoms only).'''
    from chimerax.isolde.atomic.template_utils import find_incorrect_residues
    questionable = find_incorrect_residues(session, model, heavy_atoms_only=True)
    labels = sorted(
        '{} {}{}'.format(r.name, r.chain_id, r.number) for r in questionable)
    golden('incorrect_residues', {
        'count': len(questionable),
        'labels': labels,
    })


def test_residue_graph_matches_connectivity(session, model):
    '''The perceived graph has one node per atom and one edge per intra-residue
    bond -- the substrate any graph-ML parameteriser would consume.'''
    from chimerax.isolde.atomic.template_utils import residue_graph
    # Pick a sizeable residue so the graph is non-trivial.
    r = max(model.residues, key=lambda r: len(r.atoms))
    g = residue_graph(r, label='element')

    # mcsplit.Graph: one node per atom, labels = element numbers, symmetric
    # adjacency for the intra-residue bonds.
    assert g.n == len(r.atoms)
    assert numpy.array_equal(numpy.asarray(g.labels), r.atoms.elements.numbers)
    n_edges = int(numpy.asarray(g.adjacency_matrix).sum() // 2)
    assert n_edges == len(r.atoms.intra_bonds)
