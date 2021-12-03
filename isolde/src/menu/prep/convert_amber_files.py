# @Author: Tristan Croll <tic20>
# @Date:   15-Jul-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 03-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

tooltip = ('Convert an AMBER .frcmod/.mol2 pair to the OpenMM ffXML format')

def run_script(session):
    from chimerax.core.errors import UserError
    from chimerax.isolde.parmed import install_parmed_if_necessary
    install_parmed_if_necessary(session)
    from Qt.QtWidgets import QFileDialog
    import os
    filter = "AMBER Parameter files (*.mol2 *.frcmod)"
    files, _ = QFileDialog.getOpenFileNames(None, "Choose one .mol2 and one .frcmod file", "", filter)
    if len(files)!=2:
        raise UserError("Please select exactly one .mol2 and one .frcmod file!")
    matched_files = {}
    for f in files:
        matched_files[os.path.splitext(f)[1].lower()] = f
    if '.mol2' not in matched_files.keys() or '.frcmod' not in matched_files.keys():
        raise UserError('Please select exactly one .mol2 and one .frcmod file!')
    mol2, frcmod = matched_files['.mol2'], matched_files['.frcmod']
    from chimerax.isolde.openmm.amberff.amber_convert import amber_to_ffxml
    try:
        ffxml = amber_to_ffxml(frcmod, mol2)
    except Exception as e:
        raise UserError(f'ParmEd failed to run with the following error. Do your .mol2 and .frcmod match?\n {str(e)}')
    session.logger.info(f'Converted AMBER files {mol2} and {frcmod} to OpenMM FFXML file {ffxml}.')
    if hasattr(session, 'isolde'):
        ff_name = session.isolde.sim_params.forcefield
        forcefield = session.isolde.forcefield_mgr[ff_name]
        forcefield.loadFile(ffxml, resname_prefix='USER_')
        session.logger.info('This has been loaded into ISOLDE\'s forcefield for this session. For '
            'future sessions you should add it using the "Load residue MD definition(s)" button.')
    else:
        session.logger.info('On starting ISOLDE you can add it to the forcefield using the '
            '"Load residue MD definition(s)" button.')
        