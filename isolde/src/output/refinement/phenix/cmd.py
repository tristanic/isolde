from chimerax.core.commands import (
    register, CmdDesc,
    FloatArg, IntArg, BoolArg, FileNameArg,
    EnumOf
    
)
from chimerax.atomic import AtomicStructuresArg, StructureArg
from chimerax.map import MapArg
from chimerax.core.errors import UserError

def write_rsr_cmd(session, model, resolution, volume, 
        model_file_name=None, param_file_name=None, 
        restrain_coordination_sites=False,
        restrain_positions=False,
        include_hydrogens=False,
        run_dssp=True):
    from .real_space_refine_input import write_real_space_refine_defaults
    write_real_space_refine_defaults(
        session, model, volume, resolution, 
        model_file_name=model_file_name, param_file_name=param_file_name,
        restrain_coordination_sites=restrain_coordination_sites,
        restrain_positions=restrain_positions,
        include_hydrogens=include_hydrogens,
        run_dssp=run_dssp
    )

def write_phenix_refine_cmd(session, model, 
        model_file_name=None, param_file_name=None,
        restrain_coordination_sites=False,
        include_hydrogens=False,
        num_processors=1,
        num_macrocycles=6,
        nqh_flips=True,
        scattering_type='xray',
        run_dssp=True):
    from chimerax.clipper import get_symmetry_handler
    sh = get_symmetry_handler(model, create=False)
    crystal_error_str = 'Model must be a crystal structure initialised in Clipper with experimental data!'
    if sh is None:
        raise UserError(crystal_error_str)
    if not len(sh.map_mgr.xmapsets):
        raise UserError(crystal_error_str)
    xmapset = sh.map_mgr.xmapsets[0]
    from .refine_input import write_phenix_refine_defaults
    write_phenix_refine_defaults(session, model, xmapset, 
        model_file_name=model_file_name, param_file_name=param_file_name,
        restrain_coordination_sites=restrain_coordination_sites,
        include_hydrogens=include_hydrogens,
        num_processors=num_processors,
        num_macrocycles=num_macrocycles,
        nqh_flips=nqh_flips,
        scattering_type=scattering_type,
        run_dssp=run_dssp)

def register_write_phenix_rsr(logger):
    desc = CmdDesc(
        required=[
            ('model', StructureArg),
            ('resolution', FloatArg),
            ('volume', MapArg),
        ],
        keyword=[
            ('model_file_name', FileNameArg),
            ('param_file_name', FileNameArg),
            ('restrain_positions', BoolArg),
            ('include_hydrogens', BoolArg),
            ('run_dssp', BoolArg)
        ],
        synopsis='Write input files for phenix.real_space_refine, restraining the model to its current geometry.'
    )
    register('isolde write phenixRsrInput', desc, write_rsr_cmd, logger=logger)

def register_write_phenix_refine(logger):
    desc = CmdDesc(
        required=[
            ('model', StructureArg),
        ],
        keyword=[
            ('model_file_name', FileNameArg),
            ('param_file_name', FileNameArg),
            ('include_hydrogens', BoolArg),
            ('num_processors', IntArg),
            ('num_macrocycles', IntArg),
            ('nqh_flips', BoolArg),
            ('scattering_type', EnumOf(('xray','electron','neutron'))),
            ('run_dssp', BoolArg)
        ],
        synopsis='Write input files for phenix.refine, restraining the model to its current geometry.'
    )
    register('isolde write phenixRefineInput', desc, write_phenix_refine_cmd, logger=logger)

def register_phenix_commands(logger):
    register_write_phenix_rsr(logger)
    register_write_phenix_refine(logger)