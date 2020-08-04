# @Author: Tristan Croll <tic20>
# @Date:   01-Aug-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 01-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

global _have_parmed
_have_parmed = False

def have_parmed():
    global _have_parmed
    if _have_parmed:
        return _have_parmed
    try:
        import parmed as pmd
        _have_parmed = True
        return True
    except ImportError:
        return False

def install_parmed_if_necessary(session):
    if not have_parmed():
        from chimerax.isolde.dialog import choice_warning
        msg = ("This command requires ParmEd, which needs a system compiler for "
            "installation. Would you like to try to install it now?")
        result = choice_warning(msg)
        if result:
            install_parmed(session)

def install_parmed(session):
    import sys, os, glob, subprocess
    from chimerax import app_bin_dir, app_dirs
    install_dir = app_dirs.user_data_dir
    python_executable = glob.glob(os.path.join(app_bin_dir, "python*"))[0]
    cmd_args = [python_executable, "-m", "pip", "install",
        "-t", os.path.join(app_dirs.user_data_dir,'site-packages'),
        "--upgrade", "--force-reinstall", "versioneer", "parmed"]
    pipes = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = pipes.communicate()

    if pipes.returncode != 0:
        platform = sys.platform.lower()
        if platform=='darwin':
            compiler_package = 'XCode'
        elif platform=='windows':
            compiler_package = 'Visual Studio'
        from chimerax.core.errors import UserError
        err_str = ('Installing ParmEd failed with the following error. '
            'The most common cause of failure is lack of a suitable compiler. '
            'Please install {} and try again. \n\n{}'.format(compiler_package,
            std_err.strip()))
        from chimerax.core.errors import UserError
        raise UserError(err_str)

    session.logger.info(std_out.strip().decode('utf-8'))
    session.logger.warning(std_err.strip().decode('utf-8'))
