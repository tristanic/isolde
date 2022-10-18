# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 31-Jul-2020
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
        'free_amino_acids.xml',      # Free amino acids, DOI: 10.1007/s00894-014-2478-z
        ]],

    'charmm36': ['charmm36.xml', 'charmm36/water.xml',]
}

_ligand_files = {
    'amber14': os.path.join(ff_dir, 'amberff', 'moriarty_and_case.zip'),
    'charmm36': None
}

default_forcefields = list(_forcefield_files.keys())

def _define_forcefield(ff_files):
    #from openmm.app import ForceField
    ff = ForceField(*[f for f in ff_files if f is not None])
    return ff

def _background_load_ff(name, ff_files, openmm_version, isolde_version):
    cache_dir = get_forcefield_cache_dir()
    pickle_file = os.path.join(cache_dir,'ff_{}.pickle'.format(name))
    try:
        import pickle
        with open(pickle_file, 'rb') as infile:
            forcefield, cached_openmm_version, cached_isolde_version = pickle.load(infile)
        if cached_openmm_version != openmm_version or cached_isolde_version != isolde_version:
            raise RuntimeError('Cached forcefield is out of date.')
        return {name: forcefield}
    except:
        print('Forcefield cache not found or out of date. Regenerating from ffXML files. This is normal if running ISOLDE for the first time, or after upgrading OpenMM.')
        ff = _define_forcefield(ff_files)
        with open(pickle_file, 'wb') as outfile:
            pickle.dump((ff, openmm_version, isolde_version), outfile)
        print('Done loading forcefield')
        return {name: ff}

class ForcefieldMgr:
    def __init__(self, session):
        self.session=session
        self._ff_dict = {}
        self._ligand_dict = {}
        self._task = None
        from openmm import version
        self._openmm_version = version.version
        from chimerax.isolde import __version__
        self._isolde_version = __version__

    def reset(self, clear_cache=False):
        if clear_cache:
            self.clear_cache()
        self._ff_dict.clear()

    @staticmethod
    def clear_cache():
        ff_dir = get_forcefield_cache_dir()
        import glob, os
        for f in glob.glob(os.path.join(ff_dir, '*.pickle')):
            os.remove(f)



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
                    namelist.append(name)
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

from openmm.app import ForceField as _ForceField
class ForceField(_ForceField):
    def assignTemplates(self, topology, ignoreExternalBonds=False,
            explicit_templates={}):
        '''
        Parameters
        ----------

        topology: openmm `Topology`
            The Topology whose residues are to be checked against the forcefield
            residue template
        ignoreExternalBonds: bool=False
            If true, ignore external bonds when matching residues to templates.
        explicit_templates: dict={}
            An optional {residue: template_name} dict specifying the templates to
            use for particular residues


        Returns three items:

            - a {residue: template} dict for all residues in the topology for
              which a unique template was found
            - a {residue: [templates]} dict for any residues with multiple
              matching templates
            - a list of residues with no matching template.

        For a simulation to start the second and third objects should be empty.
        '''
        from openmm.app.forcefield import _createResidueSignature
        from openmm.app.internal import compiled
        bondedToAtom = self._buildBondedToAtomList(topology)
        unique_matches = {}
        multiple_matches = {}
        unmatched = []

        templateSignatures = self._templateSignatures
        for res in topology.residues():
            sig = _createResidueSignature([atom.element for atom in res.atoms()])
            explicit = explicit_templates.get(res, None)
            if explicit:
                t = self._templates[explicit]
                match = compiled.matchResidueToTemplate(res, t, bondedToAtom, ignoreExternalBonds)
                if match is not None:
                    unique_matches[res] = (t, match)
                else:
                    unmatched.append(res)
                continue
            if sig in templateSignatures:
                allMatches = []
                for t in templateSignatures[sig]:
                    match = compiled.matchResidueToTemplate(res, t, bondedToAtom, ignoreExternalBonds)
                    if match is not None:
                        allMatches.append((t, match))
                if len(allMatches) == 1:
                    unique_matches[res] = allMatches[0]
                elif len(allMatches) > 1:
                    multiple_matches[res] = allMatches
                else:
                    unmatched.append(res)
            else:
                unmatched.append(res)
        return unique_matches, multiple_matches, unmatched

    def findNearestTemplates(self, residue):
        '''
        Find potential templates for the given residue by name and topology.
        '''
        pass

    def registerResidueTemplate(self, template):
        super().registerResidueTemplate(template)
        if len(template.atoms) > 2:
            template.graph = self.template_graph(template)
        else:
            template.graph = None

    @staticmethod
    def residue_graph(residue):
        '''
        Make a graph representing the connectivity of a residue's atoms.
        '''
        from chimerax.isolde.graph import Graph
        import numpy
        atoms = [a for a in residue.atoms()]
        labels = [a.element.atomic_number for a in atoms]
        edges = numpy.array([[atoms.index(b.atom1), atoms.index(b.atom2)] for b in residue.internal_bonds()])
        return Graph(labels, edges)

    @staticmethod
    def template_graph(template):
        '''
        Make a graph representing the connectivity of a template's atoms.
        '''
        from chimerax.isolde.graph import Graph
        import numpy
        atoms = template.atoms
        labels = [a.element.atomic_number for a in atoms]
        edges = numpy.array(list(template.bonds))
        try:
            return Graph(labels, edges)
        except:
            raise RuntimeError('Failed to make graph for {}. Labels: {}\nEdges: {}\n'.format(
                template.name, labels, edges
            ))

    @staticmethod
    def match_score(residue, template, residue_indices):
        ratoms = list(residue.atoms())
        residue_size = sum(a.element.atomic_number for a in ratoms)
        template_size = sum(a.element.atomic_number for a in template.atoms)
        match_size = sum(ratoms[i].element.atomic_number for i in residue_indices)
        residue_delta = residue_size - match_size
        template_delta = template_size - match_size
        return 1 - (residue_delta + template_delta)/min(residue_size, template_size)


    def find_possible_templates(self, res, max_missing_heavy_atoms=3,
            max_missing_heavy_atom_fraction = 0.2,
            minimum_graph_match = 0.9):
        from openmm.app import element
        from collections import Counter
        from ..atomic.template_utils import find_maximal_isomorphous_fragment
        from .amberff.glycam import (
            known_sugars,
            residue_name_to_glycam_code,
            _glycam_prefix
            )
        residue_counts = Counter(a.element for a in res.atoms())
        residue_elements = set(residue_counts.keys())
        num_atoms = sum(residue_counts.values())
        num_heavy_atoms = sum(count for e, count in residue_counts.items() if e not in (None, element.hydrogen))
        matches_by_name = []
        matches_by_composition = []
        rname = res.name
        if num_atoms >= 3:
            rgraph = self.residue_graph(res)
        else:
            rgraph = None
        if rname in known_sugars:
            base_name = residue_name_to_glycam_code[rname]
            o_links = []
            for b in res.external_bonds():
                a0, a1 = b
                if a0.residue == res:
                    if a0.name.startswith('O'):
                        o_links.append(int(a0.name[1]))
                elif a1.name.startswith('O'):
                    o_links.append(int(a1.name[1]))
            o_links = tuple(sorted(o_links))
            if not len(o_links):
                o_links = (0,)
            prefix = _glycam_prefix.get(o_links, None)

            if prefix is not None:
                tmpl_name = 'GLYCAM_'+prefix+base_name
                template = self._templates[tmpl_name]
                tgraph = template.graph
                if rgraph:
                    r_indices, t_indices, timed_out = rgraph.maximum_common_subgraph(tgraph, timeout=0.25)
                    score = self.match_score(res, template, r_indices)
                else:
                    score = -1
                matches_by_name.append((tmpl_name, score))
        for template_name, template in self._templates.items():
            split_name = template_name.split('_')
            if len(split_name) > 1:
                base_name = split_name[1].upper()
            else:
                base_name = template_name.upper()
            tgraph = template.graph
            # if len(template.atoms) >= 3:
            #     tgraph = self.template_graph(template)
            # else:
            #     tgraph = None
            # Keep all templates with similar names
            if rname.upper() == base_name:
                if rgraph and tgraph:
                    r_indices, t_indices, aborted = rgraph.maximum_common_subgraph(tgraph, timeout=0.25)
                    score = self.match_score(res, template, r_indices)
                else:
                    score = -1
                matches_by_name.append((template_name, score))
                continue
            # If the residue has at least 3 atoms, try template matching
            if rgraph is None or tgraph is None:
                continue
            template_counts = Counter(a.element for a in template.atoms)
            all_elements = set(template_counts.keys()).union(residue_elements)
            differences = {e:abs(template_counts.get(e, 0) - residue_counts.get(e, 0)) for e in all_elements}
            heavy_differences = {e: d for e, d in differences.items() if e != element.hydrogen}
            hd = sum(heavy_differences.values())
            if hd > max_missing_heavy_atoms and hd/num_heavy_atoms > max_missing_heavy_atom_fraction:
                continue
            template_size = len(template.atoms)
            r_indices, t_indices, aborted = rgraph.maximum_common_subgraph(tgraph, timeout=0.25)
            num_matched = len(r_indices)
            score = self.match_score(res, template, r_indices)
            if num_matched/num_atoms > minimum_graph_match:
                matches_by_composition.append((template_name, score))
        matches_by_name = list(sorted(matches_by_name, reverse=True, key=lambda t: t[1]))
        # matches_by_name = [t[0] for t in matches_by_name]
        matches_by_composition = list(sorted(matches_by_composition, reverse=True, key=lambda t: t[1]))
        # matches_by_name = [t[0] for t in matches_by_composition]

        return matches_by_name, matches_by_composition
