# @Author: Tristan Croll <tic20>
# @Date:   03-Aug-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 03-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

tooltip = ('Convert an ASP or GLU sidechain to the acidic form')

def run_script(session):
    from chimerax.atomic import selected_residues
    from chimerax.build_structure import modify_atom
    from chimerax.build_structure.mod import ParamError
    from chimerax.core.errors import UserError

    sel = selected_residues(session)
    if len(sel) != 1 or sel[0].name not in ('ASP', 'GLU'):
        raise UserError('Please select a single ASP or GLU residue!')
    sel = sel[0]
    if sel.name=='ASP':
        pos = 'D'
    else:
        pos = 'E'
    o_atom = sel.find_atom(f'O{pos}2')
    other_o = sel.find_atom(f'O{pos}1')
    if o_atom is None or other_o is None:
        raise UserError('Selected acid sidechain is missing one or both of its oxygen atoms!')
    if len(o_atom.neighbors) != 1 or len(other_o.neighbors) != 1:
        raise UserError('Selected acid sidechain already has a substituent on its carboxyl group!')
    
    new_h = modify_atom(o_atom, o_atom.element, 2, connect_back=False, res_name=sel.name)[1]
    new_h.color = [255,255,255,255]
