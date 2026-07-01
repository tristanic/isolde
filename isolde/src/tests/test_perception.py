# @Author: Tristan Croll
# @Date:   30-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 30-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Characterisation test for ISOLDE's existing chemical-perception helpers --
comparing modelled residues against templates, and building a connectivity graph
from a residue. The future ML parameterisation layer leans on perception of
exactly this kind. Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_perception.py

Asserts the CURRENT behaviour for the bundled 1pmx fixture. Prints PASS/FAIL and
exits non-zero on failure.
'''
import os
import numpy

FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '1pmx_1.pdb')

# Baseline captured 2026-06: every residue in 1pmx matches its template.
EXPECTED_INCORRECT_COUNT = 0


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def run(session):
    from chimerax.core.commands import run as run_cmd
    m = run_cmd(session, 'open "%s"' % FIXTURE.replace('\\', '/'), log=False)[0]
    try:
        # --- template comparison ------------------------------------------
        from chimerax.isolde.atomic.template_utils import (
            find_incorrect_residues, residue_graph)
        questionable = find_incorrect_residues(session, m, heavy_atoms_only=True)
        if len(questionable) != EXPECTED_INCORRECT_COUNT:
            labels = ['%s %s%d' % (r.name, r.chain_id, r.number) for r in questionable]
            _fail('find_incorrect_residues changed: expected %d, got %d %r'
                  % (EXPECTED_INCORRECT_COUNT, len(questionable), labels))
        print('PASS: find_incorrect_residues flags %d residues (baseline)'
              % len(questionable))

        # --- perceived connectivity graph ---------------------------------
        # mcsplit.Graph: one node per atom, labels = element numbers, symmetric
        # adjacency for the intra-residue bonds.
        r = max(m.residues, key=lambda r: len(r.atoms))
        g = residue_graph(r, label='element')
        if g.n != len(r.atoms):
            _fail('residue_graph node count %d != atom count %d' % (g.n, len(r.atoms)))
        if not numpy.array_equal(numpy.asarray(g.labels), r.atoms.elements.numbers):
            _fail('residue_graph labels do not match atom element numbers')
        n_edges = int(numpy.asarray(g.adjacency_matrix).sum() // 2)
        if n_edges != len(r.atoms.intra_bonds):
            _fail('residue_graph edge count %d != intra-bond count %d'
                  % (n_edges, len(r.atoms.intra_bonds)))
        print('PASS: residue_graph(%s) faithfully encodes %d atoms / %d bonds'
              % (r.name, g.n, n_edges))
    finally:
        session.models.close([m])

    print('ALL PASS')


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
