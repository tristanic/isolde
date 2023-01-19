# @Author: Tristan Croll <tic20>
# @Date:   17-Jul-2020
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 17-Jul-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def reset_forcefield(session):
    from chimerax.isolde.openmm.forcefields import ForcefieldMgr
    ForcefieldMgr.clear_cache()
    if hasattr(session, 'isolde'):
        session.isolde.forcefield_mgr.reset()



def register_ff_cmd(logger):
    from chimerax.core.commands import register, CmdDesc
    register('isolde reset forcefield',
        CmdDesc(synopsis='Reset the forcefield and regenerate from ffXML files'),
        reset_forcefield,
        logger=logger
    )
