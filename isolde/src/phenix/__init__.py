# @Author: Tristan Croll <tic20>
# @Date:   12-Dec-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 04-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def check_for_phenix():
    import os
    phenix_dir = os.environ.get('PHENIX', None)
    return (phenix_dir is not None)
