# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



import os

cwd = os.path.dirname(os.path.abspath(__file__))
_forcefield_files = {
    'amber14':  [os.path.join(cwd, 'amberff', f) for f in
        ['amberff14SB.xml',     # Standard AMBER14 protein/nucleic acid
        'tip3p_standard.xml',   # TIP3P water
        'tip3p_HFE_multivalent.xml', # Metal ions
        'tip3p_IOD_multivalent.xml', # Metal ions
        'gaff2.xml',                 # General AMBER force field
        # 'all_modrna08.xml',          # Naturally occurring modified RNA bases. DOI: 10.1021/ct600329w
        'glycam_all.xml',            # GLYCAM06 force field
        'combined_ccd.xml',          # General ligands (Nigel Moriarty / Dave Case)
        'ATP.xml',
        'ptms.xml',                  # Post-translational modifications from ares.tamu.ed/FFPTM DOI: 10.1021/ct400556v
        'truncated_aa.xml',          # Artifical amino acid "stubs" to support common truncations used in model building
        ]],

    'charmm36': ['charmm36.xml', 'charmm36/water.xml',]
}

default_forcefields = list(_forcefield_files.keys())

def _define_forcefield(ff_files):
    from simtk.openmm.app import ForceField
    ff = ForceField(*[f for f in ff_files if f is not None])
    return ff

class Forcefield_Mgr:
    def __init__(self):
        self._ff_dict = {}
        self._task = None

    def _complete_task(self):
        from time import sleep
        if self._task:
            while not self._task.done():
                sleep(0.01)
        self._task = None

    def __getitem__(self, key):
        self._complete_task()
        ffd = self._ff_dict
        if key in ffd.keys():
            return ffd[key]
        else:
            try:
                ff_files = _forcefield_files[key]
            except KeyError:
                raise KeyError('No forcefield with that name has been defined! '
                    'Known forcefields are: {}'.format(
                        ', '.join(set(ffd.keys()).union(_forcefield_files.keys()))
                    ))
            ffd[key] = ff = _define_forcefield(ff_files)
            return ff

    @property
    def available_forcefields(self):
        return list(set(self._ff_dict.keys()).union(_forcefield_files.keys()))

    def define_custom_forcefield(self, name, ff_files):
        '''
        Define a custom forcefield from a set of OpenMM ffXML files.
        '''
        self._complete_task()
        self._ff_dict[name] = _define_forcefield(ff_files)


    def _background_load_ff(self, name, ff_files):
        ff = _define_forcefield(ff_files)
        self._ff_dict[name] = ff

    def background_load_ff(self, name, ff_files = None):
        '''
        Prepare a forcefield in a separate thread to reduce disruption to the
        gui.
        '''
        self._complete_task()
        if ff_files is None:
            ff_files = _forcefield_files[name]
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=1) as executor:
            self._task = executor.submit(self._background_load_ff, name, ff_files)
