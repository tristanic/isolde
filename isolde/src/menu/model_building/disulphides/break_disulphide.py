# @Author: Tristan Croll <tic20>
# @Date:   01-Aug-2020
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 03-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

tooltip = ('Break a disulphide bond to make two free cysteine residues.')

def run_script(session):
    from chimerax.atomic import selected_residues, selected_bonds, Residues
    sel = selected_residues(session)
    if len(sel) != 2 or not all(sel.names=='CYS'):
        b = selected_bonds(session)
        if len(b) == 1:
            b = b[0]
            sel = Residues([a.residue for a in b.atoms])
    if len(sel) != 2 or not all(sel.names=='CYS'):
        from chimerax.core.errors import UserError
        raise UserError('Please select exactly two cysteine residues!')

    from chimerax.isolde.atomic.building.build_utils import break_disulfide
    break_disulfide(*sel)
    session.logger.info('Broke disulphide bond between {}.'.format(
    ' and '.join(['{}{}{}'.format(c.chain_id, c.number, c.insertion_code) for c in sel])
    ))
