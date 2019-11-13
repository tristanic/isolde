# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 20-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



import os

def get_forcefield_cache_dir():
    from chimerax import app_dirs
    user_data_dir = os.path.join(app_dirs.user_data_dir, 'isolde_data', 'openmm')
    if not os.path.exists(user_data_dir):
        os.makedirs(user_data_dir)
    return user_data_dir


ff_dir = os.path.dirname(os.path.abspath(__file__))
_forcefield_files = {
    'amber14':  [os.path.join(ff_dir, 'amberff', f) for f in
        ['amberff14SB.xml',     # Standard AMBER14 protein/nucleic acid
        'tip3p_standard.xml',   # TIP3P water
        'protonation_states.xml',    # Residues with unusual alternate protonation states
        'tip3p_HFE_multivalent.xml', # Metal ions
        # 'tip3p_IOD_multivalent.xml', # Metal ions
        'gaff2.xml',                 # General AMBER force field
        'CDL.xml',                   # Cardiolipin
        'lipid17.xml',               # common lipids
        # 'all_modrna08.xml',          # <BUGGY/BROKEN> Naturally occurring modified RNA bases. DOI: 10.1021/ct600329w
        'mse.xml',                   # Approximation (same charges as MET)
        'termods.xml',                   # Various chain-terminal residue modifications
        'glycam_all.xml',            # GLYCAM06 force field
        'iron_sulfur.xml',
        #'combined_ccd.xml',          # General ligands (Nigel Moriarty / Dave Case)
        #'moriarty_and_case.xml',     # General ligands (Nigel Moriarty / Dave Case)
        #'FES.xml',                   # based on SF4 from Moriarty/Case set
        'bryce_set.xml',             # A small collection of ligands and PTMs (most notably ATP/GDP/NAD(P)(H)/FAD(H)) from http://research.bmh.manchester.ac.uk/bryce/amber
        'ptms.xml',                  # Post-translational modifications from ares.tamu.ed/FFPTM DOI: 10.1021/ct400556v
        'truncated_aa.xml',          # Artifical amino acid "stubs" to support common truncations used in model building
        ]],

    'charmm36': ['charmm36.xml', 'charmm36/water.xml',]
}

_ligand_files = {
    'amber14': os.path.join(ff_dir, 'amberff', 'moriarty_and_case.zip'),
    'charmm36': None
}

default_forcefields = list(_forcefield_files.keys())

def _define_forcefield(ff_files):
    from simtk.openmm.app import ForceField
    ff = ForceField(*[f for f in ff_files if f is not None])
    return ff

def _background_load_ff(name, ff_files, openmm_version, isolde_version):
    cache_dir = get_forcefield_cache_dir()
    pickle_file = os.path.join(cache_dir,'ff_{}.pickle'.format(name))
    try:
        import pickle
        infile = open(pickle_file, 'rb')
        forcefield, cached_openmm_version, cached_isolde_version = pickle.load(infile)
        if cached_openmm_version != openmm_version or cached_isolde_version != isolde_version:
            raise RuntimeError('Cached forcefield is out of date.')
        print('Done loading forcefield')
        return {name: forcefield}
    except:
        print('Forcefield cache not found or out of date. Regenerating from ffXML files. This is normal if running ISOLDE for the first time, or after upgrading OpenMM.')
        ff = _define_forcefield(ff_files)
        outfile = open(pickle_file, 'wb')
        pickle.dump((ff, openmm_version, isolde_version), outfile)
        outfile.close()
        print('Done loading forcefield')
        return {name: ff}

class Forcefield_Mgr:
    def __init__(self, session):
        self.session=session
        self._ff_dict = {}
        self._ligand_dict = {}
        self._task = None
        from simtk.openmm import version
        self._openmm_version = version.version
        from chimerax.isolde import version
        self._isolde_version = version.version(session)

    def _complete_task(self):
        from time import sleep
        if self._task:
            while not self._task.done():
                sleep(0.01)
            result = self._task.result()
            if isinstance(result, dict):
                self._ff_dict.update(result)
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

    def ligand_db(self, key):
        db = self._ligand_dict.get(key, None)
        if db is None:
            ligand_zip = _ligand_files[key]
            if ligand_zip is not None:
                db = (ligand_zip, self._ligand_db_from_zip(ligand_zip))
                self._ligand_dict[key] = db
        return db

    def _ligand_db_from_zip(self, ligand_zip):
        from zipfile import ZipFile
        import os
        namelist = []
        with ZipFile(ligand_zip) as zf:
            for fname in zf.namelist():
                name, ext = os.path.splitext(fname)
                if ext.lower() == '.xml':
                    namelist.append(fname)
        return namelist




    @property
    def available_forcefields(self):
        self._complete_task()
        return list(set(self._ff_dict.keys()).union(_forcefield_files.keys()))

    @property
    def loaded_forcefields(self):
        self._complete_task()
        return list(self._ff_dict.keys())

    def define_custom_forcefield(self, name, ff_files):
        '''
        Define a custom forcefield from a set of OpenMM ffXML files.
        '''
        self._complete_task()
        self._ff_dict[name] = _define_forcefield(ff_files)

    def background_load_ff(self, name, ff_files = None):
        '''
        Prepare a forcefield in a separate thread to reduce disruption to the
        gui.
        '''
        self._complete_task()
        if ff_files is None:
            ff_files = _forcefield_files[name]
        from concurrent.futures import ThreadPoolExecutor
        executor = ThreadPoolExecutor(max_workers=1)
        self._task = executor.submit(_background_load_ff, name, ff_files, self._openmm_version, self._isolde_version)
        executor.shutdown(wait=False)
