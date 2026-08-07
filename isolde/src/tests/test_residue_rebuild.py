# @Author: Tristan Croll
# @Date:   15-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 15-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Headless regression test for the residue-rebuild geometry contract
(``chimerax.isolde.atomic.template_utils.fix_residue_from_template``):

    Rebuilding a residue from its coordinate template must NOT move existing
    heavy-atom geometry unless a genuinely new heavy atom needs accommodating.

The rebuild always strips and re-adds hydrogens. A past regression let those
rebuilt hydrogens count as "freshly built", so the pendant-torsion optimiser saw
a built atom on both sides of nearly every rotatable bond (every heavy atom
carries an H) and swung intact heavy-atom arms -- e.g. a heme propionate -- into
a new pose to chase density/clashes. The fix restricts the optimiser's trigger to
freshly-built *heavy* atoms.

This test rebuilds glycerol from heavy atoms only, WITH a deliberate steric clash
against one arm (the exact condition that drove the old rotation), and asserts the
heavy skeleton does not move while the hydrogens are correctly re-added. Run:

    run_chimerax.bat --nogui --exit --script src/tests/test_residue_rebuild.py
'''


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def run(session):
    import numpy
    from chimerax import mmcif
    from chimerax.atomic import AtomicStructure
    from chimerax.atomic.struct_edit import add_bond
    from chimerax.isolde.atomic.template_utils import fix_residue_from_template

    # A small ligand with rotatable arms whose heavy atoms all carry H -- the case
    # that used to trip the optimiser. Fall back through a few common CCD ids in
    # case one is unavailable in the local component dictionary.
    tmpl = None
    for ccd_id in ('GOL', 'EDO', 'PEG'):
        try:
            tmpl = mmcif.find_template_residue(session, ccd_id)
        except Exception:
            tmpl = None
        if tmpl is not None:
            break
    if tmpl is None:
        print('SKIP: no suitable CCD template available')
        return

    s = AtomicStructure(session, name='rebuild-test', auto_style=False)
    try:
        r = s.new_residue(tmpl.name, 'A', 1)
        made = {}
        for ta in tmpl.atoms:
            if ta.element.number == 1:
                continue  # heavy only; the rebuild must re-add H without moving these
            a = s.new_atom(ta.name, ta.element)
            a.coord = numpy.array(ta.coord)
            r.add_atom(a)
            made[ta.name] = a
        for ta in tmpl.atoms:
            a1 = made.get(ta.name)
            if a1 is None:
                continue
            for nb in ta.neighbors:
                a2 = made.get(nb.name)
                if a2 is not None and a2 not in a1.neighbors:
                    add_bond(a1, a2)
        if len(made) < 3:
            print('SKIP: template %s has too few heavy atoms' % tmpl.name)
            return

        # Deliberate clash ~1 A from the most distal heavy atom: the old code would
        # swing that arm to relieve it; the fixed code must not touch heavy geometry.
        term = max(made.values(), key=lambda a: float(numpy.linalg.norm(a.coord)))
        clashr = s.new_residue('CLS', 'B', 1)
        ca = s.new_atom('CL', 'Cl')
        ca.coord = numpy.array(term.coord) + numpy.array([1.0, 0.0, 0.0])
        clashr.add_atom(ca)

        before = {a.name: numpy.array(a.coord) for a in r.atoms if a.element.number > 1}
        fix_residue_from_template(r, tmpl)
        after = {a.name: numpy.array(a.coord) for a in r.atoms if a.element.number > 1}

        if set(before) != set(after):
            _fail('heavy-atom set changed during rebuild: %s -> %s'
                  % (sorted(before), sorted(after)))
        maxd = max(float(numpy.linalg.norm(after[nm] - before[nm])) for nm in before)
        if maxd > 1e-3:
            _fail('rebuild moved existing heavy geometry by %.3f A (only H were '
                  'rebuilt; heavy atoms must not move)' % maxd)
        n_h = sum(1 for a in r.atoms if a.element.number == 1)
        if n_h == 0:
            _fail('rebuild did not re-add any hydrogens')
        print('PASS: rebuild preserves heavy geometry (max %.4f A) and re-adds %d H'
              % (maxd, n_h))
    finally:
        session.models.close([s])

    # --- Part 2: rebuilding a residue must NOT move a BONDED NEIGHBOUR ----------------
    # Bond.side_atoms() follows the whole structure graph across residue boundaries, so a
    # backbone/link bond's "side" is the rest of the chain. The pendant-torsion optimiser
    # must never rotate a side that leaves the residue (it would drag every downstream
    # residue -- the 1crx DNA-chain bug). Here: a glycerol whose terminal O3 is bonded to a
    # separate residue, a clash placed by that arm, and a heavy atom (C3, on the neighbour
    # side) deleted so a rebuild + torsion-opt is triggered. The neighbour must not move.
    s2 = AtomicStructure(session, name='rebuild-neighbour-test', auto_style=False)
    try:
        r = s2.new_residue(tmpl.name, 'A', 1)
        made = {}
        for ta in tmpl.atoms:
            if ta.element.number == 1:
                continue
            a = s2.new_atom(ta.name, ta.element)
            a.coord = numpy.array(ta.coord); r.add_atom(a); made[ta.name] = a
        for ta in tmpl.atoms:
            a1 = made.get(ta.name)
            if a1 is None:
                continue
            for nb in ta.neighbors:
                a2 = made.get(nb.name)
                if a2 is not None and a2 not in a1.neighbors:
                    add_bond(a1, a2)
        # Pick a terminal heavy atom (one heavy neighbour) as the link point, and the
        # heavy atom bonded to it as the one to delete+rebuild (so the rebuilt atom sits
        # on the side that crosses into the neighbour residue).
        link = None
        for a in made.values():
            heavy_nb = [n for n in a.neighbors if n.element.number > 1]
            if len(heavy_nb) == 1:
                link, inner = a, heavy_nb[0]
                break
        if link is None:
            print('SKIP: no terminal heavy atom to attach a neighbour')
        else:
            nbr_res = s2.new_residue('NBR', 'B', 1)
            nbr = s2.new_atom('N1', 'N')
            nbr.coord = numpy.array(link.coord) + numpy.array([1.3, 0.2, 0.0])
            nbr_res.add_atom(nbr)
            add_bond(link, nbr)                       # link into the neighbour residue
            # Clash right beside the inner atom's arm, to tempt the optimiser to rotate.
            clr = s2.new_residue('CLS', 'C', 1)
            cl = s2.new_atom('CL', 'Cl')
            cl.coord = numpy.array(inner.coord) + numpy.array([0.6, 0.6, 0.0])
            clr.add_atom(cl)
            before_nbr = numpy.array(nbr.coord)
            inner.delete()                            # rebuild will re-add it (built heavy)
            fix_residue_from_template(r, tmpl)
            after_nbr = numpy.array(nbr.coord)
            moved = float(numpy.linalg.norm(after_nbr - before_nbr))
            if moved > 1e-6:
                _fail('rebuilding a residue moved a BONDED NEIGHBOUR by %.3f A '
                      '(rebuild must never move atoms outside the residue)' % moved)
            print('PASS: rebuild leaves a bonded neighbour untouched (moved %.4f A)' % moved)
    finally:
        session.models.close([s2])

    print('ALL PASS')


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
