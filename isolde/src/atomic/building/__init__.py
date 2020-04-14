# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 14-Apr-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def set_new_atom_style(session, atoms):
    from chimerax.atomic import selected_atoms, Atom
    from chimerax.core.commands import run
    atoms.draw_modes = Atom.STICK_STYLE
    atoms.unique_residues.ribbon_hide_backbones=False
    current_sel = selected_atoms(session)
    session.selection.clear()
    atoms.selected = True
    run(session, "color sel bychain; color sel byhetero")
    session.selection.clear()
    current_sel.selected = True
