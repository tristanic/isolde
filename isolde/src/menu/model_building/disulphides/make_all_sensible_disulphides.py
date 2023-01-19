# @Author: Tristan Croll <tic20>
# @Date:   01-Aug-2020
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 03-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

tooltip = ('Create disulphide bonds between any free cysteine residues with disulphide-like geometry')

def run_script(session):
    from chimerax.core.commands import run
    from chimerax.core.errors import UserError
    from chimerax.isolde.atomic.building.build_utils import current_and_possible_disulfides, create_disulfide
    run(session, 'isolde start', log=False)
    m = session.isolde.selected_model
    if m is None:
        raise UserError('Select a model in ISOLDE first!')
    current, possible, ambiguous = current_and_possible_disulfides(m)
    for cys_pair in possible:
        create_disulfide(*cys_pair)
    if len(possible):
        session.logger.info('Created disulfide bonds between the following residues: \n{}'.format(
            '; '.join(['-'.join(['{}{}{}'.format (c.chain_id, c.number, c.insertion_code) for c in p]) for p in possible])
        ))
    if len(ambiguous):
        warn_str = ('The following cysteine residues are clustered too close to '
            'automatically assign disulphide-bonded pairs. Please check manually.\n{}').format(
                '\n'.join(', '.join(['{}{}{}'.format(c.chain_id, c.number,c.insertion_code) for c in amb_set]) for amb_set in ambiguous
            ))
        session.logger.warning(warn_str)
