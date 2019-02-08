# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



import os

cwd = os.path.dirname(os.path.abspath(__file__))
forcefields = {
    'amber14':  [os.path.join(cwd, 'amberff', f) for f in
        ['amberff14SB.xml', 'tip3p_standard.xml', 'tip3p_HFE_multivalent.xml',
        'tip3p_IOD_multivalent.xml', 'gaff2.xml', 'combined_ccd.xml']],

    'charmm36': ['charmm36.xml', 'charmm36/water.xml',]
}
