# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 05-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def copy_atom_style_from(session, atoms, ref_residue):
    from chimerax.atomic import selected_atoms, Atom
    from chimerax.core.commands import run
    r = ref_residue
    atoms.draw_modes = r.atoms[0].draw_mode
    residues = atoms.unique_residues
    residues.ribbon_displays=r.ribbon_display
    residues.ribbon_hide_backbones=r.ribbon_hide_backbone
    current_sel = selected_atoms(session)
    session.selection.clear()
    atoms.selected = True
    ref_c = r.atoms[r.atoms.element_names=='C']
    if len(ref_c):
        atoms[atoms.element_names=='C'].colors = ref_c[0].color
    else:
        run(session, "color sel bychain", log=False)
    run(session, "color sel byhetero", log=False)
    session.selection.clear()
    current_sel.selected = True

