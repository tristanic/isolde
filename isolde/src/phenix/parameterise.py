# @Author: Tristan Croll <tic20>
# @Date:   04-Aug-2020
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 17-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def check_if_params_exist(session, residue):
    from chimerax.atomic import Residues
    from chimerax.core.errors import UserError
    from chimerax.isolde.openmm.openmm_interface import (
        create_openmm_topology,
        find_residue_templates,
    )
    ffmgr = session.isolde.forcefield_mgr
    ff_name = session.isolde.sim_params.forcefield
    ff = ffmgr[ff_name]
    #ligand_db = ffmgr.ligand_db(ff_name)
    template_dict = find_residue_templates(Residues([residue]), ff, logger=session.logger)
    err_str_base = ('Residue name {} already maps to MD template {}. If this is a '
        'different chemical species, you will need to choose a new name.')
    if len(template_dict):
        err_str = err_str_base.format(residue.name, template_dict[0]
            )
        raise UserError(err_str)
    top, residue_templates = create_openmm_topology(residue.atoms, {})
    for r in top.residues():
        break
    assigned, ambiguous, unmatched = ff.assignTemplates(top,
        ignoreExternalBonds=False)
    if len(assigned):
        err_str = err_str_base.format(residue.name, assigned[r][0].name)
        raise UserError(err_str)

def parameterise_residue_with_elbow(residue, include_hydrogens=True, keep_elbow_files=False):
    session = residue.session
    resname = residue.name
    from chimerax.core.errors import UserError
    from chimerax import mmcif
    from . import check_for_phenix
    from chimerax.isolde.parmed import install_parmed_if_necessary
    install_parmed_if_necessary(session)
    if not check_for_phenix(session):
        raise UserError('Automatic MD parameterisaton requires an up-to-date version of the Phenix suite (https://www.phenix-online.org)')
    if len(residue.neighbors):
        raise UserError('Automatic MD parameterisation is currently only supported for free, non-covalently bound ligands.')
    if any(residue.atoms.elements.is_metal):
        raise UserError('Automatic MD parameterisation is not yet supported for metal complexes.')

    from ..cmd import isolde_start
    isolde_start(session)
    check_if_params_exist(session, residue)

    import os
    cwd = os.getcwd()
    from ..util import SafeTempDir
    with SafeTempDir():
        session.logger.status('Running eLBOW. This may take a while...')
        filename = resname + '_for_elbow.pdb'
        atoms = residue.atoms
        if not include_hydrogens:
            atoms = atoms[atoms.element_names != 'H']
        from chimerax.pdb import save_pdb
        session.selection.clear()
        atoms.selected=True
        save_pdb(session, filename, selected_only=True)
        _run_elbow(session, filename)
        expected_files = (resname+'.frcmod', resname+'.mol2')
        if not all(os.path.isfile(efile) for efile in expected_files):
            raise UserError('phenix.elbow finished without an error code, but the expected files were not written. Bailing out.')
        from chimerax.isolde.openmm.amberff.amber_convert import amber_to_ffxml
        amber_to_ffxml(*expected_files)
        xmlfile = resname+'.xml'
        if not os.path.isfile(xmlfile):
            raise UserError('Failed to convert .frcmod/.mol2 to ffXML.')
        import shutil
        if keep_elbow_files:
            import glob
            all_files = glob.glob('*')
            for f in all_files:
                shutil.copy(f, os.path.join(cwd, f))
        else:
            final_loc = os.path.join(cwd, xmlfile)
            shutil.copy(xmlfile, final_loc)
    ff_name = session.isolde.sim_params.forcefield
    ff = session.isolde.forcefield_mgr[ff_name]
    session.isolde.add_ffxml_files(ff, [xmlfile])

    info_str = ('Wrote AMBER parameters for residue {} to {}. These have been '
        'loaded into your forcefield for this session. For future ISOLDE sessions, '
        'load this file using the "Load residue MD definition(s)" button on '
        'the ISOLDE panel.').format(resname, final_loc)
    session.logger.info(info_str)






def _run_elbow(session, coord_file):
    from ..settings import basic_settings
    phenix_path = basic_settings.phenix_base_path
    import os, subprocess
    cmd_args = [os.path.join(phenix_path, "phenix.elbow"), "--amber", coord_file]
    pipes = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
)
    std_out, std_err = pipes.communicate()

    if pipes.returncode !=0:
        from chimerax.core.errors import UserError
        err_str = (f"Attempt to run phenix.elbow failed with the following error: \n{std_err.strip().decode('utf-8')}\n\n"
            'This could be due to poor geometry in your input '
            'coordinates, or chemistry that eLBOW does not understand (it often fails on charged groups). At '
            'present ISOLDE cannot help with this - you may need to install and run ANTECHAMBER for yourself.')
        raise UserError(err_str)
