# @Author: Tristan Croll <tic20>
# @Date:   07-Apr-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 03-Aug-2020
# @License: Lesser GNU Public License version 3.0 (see LICENSE.md)
# @Copyright: 2016-2019 Tristan Croll

tooltip = ('Create a bond between two atoms in your model.')

def run_script(session):
    from chimerax.atomic import selected_atoms
    from chimerax.core.errors import UserError
    sel = selected_atoms(session)
    if len(sel) != 2:
        raise UserError('Must have exactly two atoms selected!')
    us = sel.unique_structures
    if len(us) != 1:
        raise UserError('Both atoms must be from the same structure!')
    m = us[0]
    m.new_bond(*sel)
