# @Author: Tristan Croll <tic20>
# @Date:   04-Aug-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 04-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def parameterise_residue_with_elbow(residue, include_hydrogens=True, keep_elbow_files=True):
    session = residue.session
    resname = residue.name
    from chimerax.core.errors import UserError
    from chimerax import mmcif
    from . import check_for_phenix
    from chimerax.isolde.parmed import install_parmed_if_necessary
    install_parmed_if_necessary(session)
    if not check_for_phenix():
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
    import tempfile
    tmpdir = tempfile.TemporaryDirectory()
    tmpdirname = tmpdir.name
    os.chdir(tmpdirname)

    from chimerax.mol2 import write_mol2
    session.logger.status('Running eLBOW. This may take a while...')
    filename = resname + '_for_elbow.mol2'
    atoms = residue.atoms
    if not include_hydrogens:
        atoms = atoms[atoms.element_names != 'H']
    try:
        write_mol2(session, filename, atoms=residue.atoms)
        _run_elbow(session, filename)
    except:
        os.chdir(cwd)
        raise
    expected_files = (resname+'.frcmod', resname+'.mol2')
    if not all(os.path.isfile(efile) for efile in expected_files):
        os.chdir(cwd)
        raise UserError('phenix.elbow finished without an error code, but the expected files were not written. Bailing out.')
    from chimerax.isolde.openmm.amberff.amber_convert import amber_to_ffxml
    amber_to_ffxml(*expected_files)
    xmlfile = resname+'.xml'
    if not os.path.isfile(xmlfile):
        os.chdir(cwd)
        raise UserError('Failed to convert .frcmod/.mol2 to ffXML.')
    import shutil
    if keep_elbow_files:
        import glob
        all_files = glob.glob('*')
        for f in all_files:
            shutil.copy(os.path.join(tmpdirname, f), os.path.join(cwd, f))
    else:
        final_loc = os.path.join(cwd, xmlfile)
        shutil.copy(os.path.join(tmpdirname, xmlfile), final_loc)
    os.chdir(cwd)
    ff_name = session.isolde.sim_params.forcefield
    ff = session.isolde.forcefield_mgr[ff_name]
    session.isolde.add_ffxml_files(ff, [xmlfile])

    info_str = ('Wrote AMBER parameters for residue {} to {}. These have been '
        'loaded into your forcefield for this session. For future ISOLDE sessions, '
        'load this file using the "Load residue MD definition(s)" button on '
        'the ISOLDE panel.').format(resname, final_loc)
    session.logger.info(info_str)

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





def _run_elbow(session, mol2_file):
    import sys, os, subprocess
    cmd_args = ["phenix.elbow", "--amber", mol2_file]
    pipes = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = pipes.communicate()

    if pipes.returncode !=0:
        from chimerax.core.errors import UserError
        err_str = "Attempt to run phenix.elbow failed with the following error: \n{}".format(std_err.strip().decode('utf-8'))
        raise UserError(err_str)
