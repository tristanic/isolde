# @Author: Tristan Croll <tic20>
# @Date:   12-Dec-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 25-Sep-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

DEFAULT_MIN_PHENIX_VERSION=3964

def check_for_phenix_version():
    from chimerax.core.errors import UserError
    import sys, os, subprocess
    cmd_args = ["phenix.version"]
    try:
        pipes = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        raise UserError('The Phenix plugin requires Phenix to be installed and on your execution path!')
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



def check_for_phenix():
    version = check_for_phenix_version()
    return version >= DEFAULT_MIN_PHENIX_VERSION
