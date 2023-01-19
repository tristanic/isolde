# @Author: Tristan Croll <tic20>
# @Date:   03-Aug-2020
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 03-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

tooltip = ('Add a hydrogen to the currently selected atom, if chemically sensible.')

def run_script(session):
    from chimerax.atomic import selected_atoms
    from chimerax.build_structure import modify_atom
    from chimerax.build_structure.mod import ParamError
    from chimerax.core.errors import UserError

    sel = selected_atoms(session)
    if len(sel) != 1:
        raise UserError('Please select a single atom!')
    sel = sel[0]
    current_num_bonds = len(sel.neighbors)
    current_color = sel.color
    try:
        modify_atom(sel, sel.element, current_num_bonds+1, connect_back=False, res_name=sel.residue.name)
        sel.color = current_color
    except ParamError as e:
        # If modify_atom throws an error at the previous step, it will have deleted
        # any attached hydrogens and not put them back. We need to put them back here.
        modify_atom(sel, sel.element, current_num_bonds, connect_back=False, res_name=sel.residue.name)
        sel.color = current_color
        raise UserError(str(e))
