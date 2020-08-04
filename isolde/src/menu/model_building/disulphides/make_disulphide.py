# @Author: Tristan Croll <tic20>
# @Date:   01-Aug-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 03-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

tooltip = ('Make a disulphide bond between two cysteine residues')

def run_script(session):
    from chimerax.atomic import selected_residues
    from chimerax.isolde.atomic.building.build_utils import create_disulfide
    sel = selected_residues(session)
    if len(sel) != 2 or not all(sel.names=='CYS'):
        from chimerax.core.errors import UserError
        raise UserError('Please select two cysteine residues!')
    if sel[0] not in sel[1].neighbors:
        create_disulfide(*sel)
        session.logger.info('Created disulphide bond between {}.'.format(
        ' and '.join(['{}{}{}'.format(c.chain_id, c.number, c.insertion_code) for c in sel])
        ))
