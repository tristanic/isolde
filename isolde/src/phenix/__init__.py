# @Author: Tristan Croll <tic20>
# @Date:   12-Dec-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 25-Sep-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

DEFAULT_MIN_PHENIX_VERSION=3964


def check_for_phenix_version(session):
    from chimerax.core.errors import UserError
    import sys, os, subprocess
    from ..settings import basic_settings
    phenix_path = basic_settings.phenix_base_path
    if phenix_path is None:
        phenix_path = _choose_phenix_directory(session)
        basic_settings.phenix_base_path = phenix_path
    cmd_args = [os.path.join(phenix_path, "phenix.version")]
    try:
        pipes = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        session.logger.warning('Phenix installation appears to have moved or been deleted. Please provide a new path.')
        basic_settings.phenix_base_path = None
        return check_for_phenix_version(session)
    std_out, std_err = pipes.communicate()
    lines = [line.strip() for line in std_out.decode('utf-8').split('\n')]
    version = _parse_version(lines)
    return version

def _parse_version(lines):
    version = -1
    for line in lines:
        if line.startswith('Release'):
            version = int(line.split(':')[-1].strip())
    return version



def check_for_phenix(session):
    version = check_for_phenix_version(session)
    return version >= DEFAULT_MIN_PHENIX_VERSION

def _choose_phenix_directory(session):
    satisfied = False
    from Qt.QtWidgets import QFileDialog, QMessageBox
    parent = session.ui.main_window
    import subprocess, os
    while not satisfied:
        result = QFileDialog.getExistingDirectory(parent, 'Please provide the directory containing the Phenix executables.', options=QFileDialog.Options())
        if not result:
            break
        try:
            subprocess.call([os.path.join(result,'phenix.version')])
            satisfied = True
        except FileNotFoundError:
            choice = QMessageBox.warning(parent, 'This directory does not appear to contain Phenix executables. Would you like to try again?',
                QMessageBox.Ok|QMessageBox.Cancel)
        except:
            raise
    if not satisfied:
        from chimerax.core.errors import UserError
        raise UserError('Could not find Phenix installation. Operation cancelled')
    return result
