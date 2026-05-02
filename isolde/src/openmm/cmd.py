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
        session.isolde.forcefield_mgr.reset(clear_cache=True)


def _residue_names_from_ffxml(xml_file, resname_prefix=''):
    '''Parse an OpenMM ffXML file and return the list of residue names it
    defines, optionally prefixed.'''
    import xml.etree.ElementTree as ET
    tree = ET.parse(xml_file)
    root = tree.getroot()
    names = []
    residues_section = root.find('Residues')
    if residues_section is not None:
        for residue in residues_section.findall('Residue'):
            if 'name' in residue.attrib:
                names.append(resname_prefix + residue.attrib['name'])
    return names


def _delete_residue_templates(forcefield, template_names):
    '''Remove residue templates from an OpenMM ForceField, including their
    entries in the topology-signature index used for fallback matching.'''
    for name in template_names:
        if name in forcefield._templates:
            template = forcefield._templates[name]
            del forcefield._templates[name]
            from openmm.app.forcefield import _createResidueSignature
            signature = _createResidueSignature([atom.element for atom in template.atoms])
            if signature in forcefield._templateSignatures:
                if template in forcefield._templateSignatures[signature]:
                    forcefield._templateSignatures[signature].remove(template)


def _load_user_ffxml(forcefield, file_path):
    '''Load an OpenMM ffXML file into ``forcefield`` with the ``USER_`` prefix,
    first removing any existing ``USER_``-prefixed templates whose names would
    otherwise collide. Returns the list of residue names defined in the file
    (without the ``USER_`` prefix).'''
    residue_names = _residue_names_from_ffxml(file_path)
    template_names = ['USER_' + n for n in residue_names]
    _delete_residue_templates(forcefield, template_names)
    forcefield.loadFile(file_path, resname_prefix='USER_')
    return residue_names


def isolde_load_parameters(session, file_path):
    '''
    Load OpenMM ffXML residue parameters into ISOLDE's active forcefield.
    Mirrors the "Load residue parameters" GUI button: any existing USER_
    templates with colliding names are replaced.
    '''
    from chimerax.isolde.cmd.cmd import isolde_start, block_if_sim_running
    from chimerax.core.errors import UserError
    import os

    block_if_sim_running(session)
    isolde = isolde_start(session)
    if isolde is None:
        raise UserError('isolde load parameters requires the ISOLDE GUI session.')

    file_path = os.path.expanduser(file_path)
    if not os.path.isfile(file_path):
        raise UserError(f'No such file: {file_path}')

    ff = isolde.forcefield_mgr[isolde.sim_params.forcefield]
    try:
        residue_names = _load_user_ffxml(ff, file_path)
    except ValueError as e:
        session.logger.warning(f'Failed to load {file_path}: {e}')
        return

    if residue_names:
        formatted = ', '.join(residue_names)
        session.logger.info(
            f'ISOLDE: loaded residue parameters from {file_path} '
            f'({len(residue_names)} residue{"s" if len(residue_names) != 1 else ""} '
            f'registered as USER_): {formatted}'
        )
    else:
        session.logger.warning(
            f'ISOLDE: loaded {file_path}, but it contained no <Residue> definitions.'
        )


def register_ff_cmd(logger):
    from chimerax.core.commands import register, CmdDesc, FileNameArg
    register('isolde reset forcefield',
        CmdDesc(synopsis='Reset the forcefield and regenerate from ffXML files'),
        reset_forcefield,
        logger=logger
    )
    register('isolde load parameters',
        CmdDesc(
            required=[('file_path', FileNameArg)],
            synopsis='Load OpenMM ffXML residue parameters into ISOLDE'
        ),
        isolde_load_parameters,
        logger=logger,
    )
