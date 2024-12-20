
class ProblemAggregator:

    def __init__(self, session):
        self.session = session
        self._restraint_problem_getters = {
            'Standard torsion restraints': self.get_torsion_restraint_problems,
            'Standard distance restraints': self.get_distance_restraint_problems,
            'Adaptive torsion restraints': self.get_adaptive_torsion_restraint_problems,
            'Reference distance restraints': self.get_adaptive_distance_restraint_problems,
        }

        self._validation_problem_getters = {
            'Protein sidechains':               self.get_rotamer_problems,
            'Protein backbone':                 self.get_protein_backbone_problems,
            'Clashes':                          self.get_clashes,
        }

    from chimerax.isolde.molobject import (
        ProperDihedralRestraint, DistanceRestraint,
        AdaptiveDihedralRestraint, AdaptiveDistanceRestraint,
        Rotamer, Rama
    )
    from chimerax.isolde.validation.clashes import Clash

    _registered_types = {
        'Standard torsion restraints': ProperDihedralRestraint,
        'Standard distance restraints': DistanceRestraint,
        'Adaptive torsion restraints': AdaptiveDihedralRestraint,
        'Reference distance restraints': AdaptiveDistanceRestraint,
        'Protein sidechains': Rotamer,
        'Protein backbone': Rama,
        'Clashes': Clash
    }

    _registered_names = {val:key for key, val in _registered_types.items()}


    @property
    def registered_restraint_problem_types(self):
        return list(self._restraint_problem_getters.keys())
    
    @property
    def registered_validation_problem_types(self):
        return list(self._validation_problem_getters.keys())
    
    def registered_type(self, name):
        return self._registered_types[name]
    
    def registered_name(self, vtype):
        return self._registered_names[vtype]
    
    def problem_zones(self, structure, restraint_types=None, validation_types=None, cutoff=3, min_points=6, validation_outliers_only=False):
        if restraint_types is None:
            restraint_types = self.registered_restraint_problem_types
        if validation_types is None:
            validation_types = self.registered_validation_problem_types
        import numpy
        sites = []
        for restraint_type in restraint_types:
            pm = self._restraint_problem_getters.get(restraint_type, None)
            if pm is None:
                self.session.logger.warning(f'{restraint_type} is not a type registered for problem aggregation!')
            else:
                sites.extend(pm(structure))
        for validation_type in validation_types:
            vm = self._validation_problem_getters.get(validation_type, None)
            if vm is None:
                self.session.logger.warning(f'{restraint_type} is not a type registered for problem aggregation!')
            else:
                sites.extend(vm(structure, outliers_only=validation_outliers_only))
        coords = numpy.array([site.center for site in sites])
        from ..clustering import dbscan
        clusters, noise = dbscan(coords, cutoff, min_points)
        from operator import itemgetter
        clustered_issues = []
        for c in clusters:
            f = itemgetter(*c)
            clustered_issues.append(f(sites))
        if len(noise):
            f = itemgetter(*noise)
            remainder = f(sites)
        else:
            remainder = []
        return clustered_issues, remainder
    
    def cluster_atoms(self, items):
        from chimerax.atomic import concatenate, Atoms
        return concatenate([Atoms(i.atoms) for i in items]).unique()

    @staticmethod
    def get_torsion_restraint_problems(structure):
        from chimerax.isolde import session_extensions as sx
        pdrm = sx.get_proper_dihedral_restraint_mgr(structure, create=True)
        torsions = pdrm.get_all_restraints_for_residues(structure.residues)
        torsions = torsions[torsions.enableds]
        return torsions[torsions.unsatisfieds]
    
    @staticmethod
    def get_distance_restraint_problems(structure):
        from chimerax.isolde import session_extensions as sx
        drm = sx.get_distance_restraint_mgr(structure, create=True)
        distances = drm.all_restraints
        distances = distances[distances.enableds]
        return distances[distances.unsatisfieds]
    
    @staticmethod
    def get_adaptive_torsion_restraint_problems(structure):
        from chimerax.isolde import session_extensions as sx
        adrm = sx.get_adaptive_dihedral_restraint_mgr(structure, create=True)
        torsions = adrm.get_all_restraints_for_residues(structure.residues)
        torsions = torsions[torsions.enableds]
        return torsions[torsions.unsatisfieds]
    
    @staticmethod
    def get_adaptive_distance_restraint_problems(structure):
        from chimerax.isolde import session_extensions as sx
        from chimerax.isolde.restraints.restraint_utils import DEFAULT_ADAPTIVE_RESTRAINT_GROUP_NAME
        drm = sx.get_adaptive_distance_restraint_mgr(structure, create=True, name=DEFAULT_ADAPTIVE_RESTRAINT_GROUP_NAME)
        distances = drm.all_restraints
        distances = distances[distances.enableds]
        return distances[distances.unsatisfieds]

    @staticmethod
    def get_rotamer_problems(structure, outliers_only=True):
        from chimerax.isolde import session_extensions as sx
        rmgr = sx.get_rotamer_mgr(structure.session)
        rotamers = rmgr.get_rotamers(structure.residues)
        if outliers_only:
            f = rmgr.outliers
        else:
            f = rmgr.non_favored_rotamers
        return f(rotamers)[0]
    
    @staticmethod
    def get_protein_backbone_problems(structure, outliers_only=True):
        from chimerax.isolde import session_extensions as sx
        rmgr = sx.get_ramachandran_mgr(structure.session)
        residues = structure.residues
        if outliers_only:
            f = rmgr.outliers
        else:
            f = rmgr.non_favored
        problems = f(residues)
        # Also count twisted and cis-nonPro peptide bonds here
        cis = rmgr.cis(residues)
        cis_nonpro = cis[cis.names != 'PRO']
        # Note: this will probably lead to double-counting since 
        # twisted peptides will usually also be captured by 
        # get_torsion_restraint_problems()... but ultimately I 
        # don't think that's a huge problem.
        twisted = rmgr.twisted(residues)
        from chimerax.atomic import concatenate, Residues
        twisted = Residues([t[0] for t in twisted])
        problems = concatenate([problems, cis_nonpro,twisted])
        problem_ramas = rmgr.get_ramas(problems)
        return problem_ramas

    @staticmethod
    def get_clashes(structure, outliers_only=True):
        from chimerax.isolde.validation.clashes import unique_clashes
        clashes = unique_clashes(structure.session, structure.atoms, severe_only=outliers_only)
        return clashes
