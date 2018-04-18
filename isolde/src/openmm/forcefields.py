# @Author: Tristan Croll
# @Date:   07-Mar-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   Tristan Croll
# @Last modified time: 18-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



import os

cwd = os.path.dirname(os.path.abspath(__file__))
forcefields = {
    'amber14':  [os.path.join(cwd, 'amberff', f) for f in
        ['amberff14SB.xml','tip3p_standard.xml',
            'tip3p_HFE_multivalent.xml', 'tip3p_IOD_multivalent.xml']],
}
