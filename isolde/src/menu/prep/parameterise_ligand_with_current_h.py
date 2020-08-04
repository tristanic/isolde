# @Author: Tristan Croll <tic20>
# @Date:   04-Aug-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 04-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

tooltip = ('Parameterise the selected residue for MD, using the current hydrogens if sensible.')

def run_script(session):
    from chimerax.atomic import selected_residues
    from chimerax.core.errors import UserError
    sel = selected_residues(session)
    if len(sel) != 1:
        raise UserError('Please select a single residue!')
    sel = sel[0]
    from chimerax.isolde.phenix.parameterise import parameterise_residue_with_elbow
    parameterise_residue_with_elbow(sel, include_hydrogens=True, keep_elbow_files=False)
