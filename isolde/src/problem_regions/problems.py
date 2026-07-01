
# Variance-normalised excess-stress (von Mises z-score) above which a live atom is
# treated as a local-strain "problem". Follows the Problem Zones "Outliers only"
# convention shared (approximately) across validation metrics: checked ~ Z=3,
# unchecked ~ Z=2. Kept fairly permissive on purpose -- the downstream DBSCAN
# min_points gate suppresses isolated flags, so the threshold need not pre-filter
# hard. get_local_strain_problems picks between these on the outliers_only flag.
STRAIN_Z_OUTLIER = 3.0      # "Outliers only" checked
STRAIN_Z_NONFAVORED = 2.0   # "Outliers only" unchecked


class StrainSite:
    '''
    Adapter presenting one live high-strain atom as a Problem-Zones "site"
    (``.center`` + ``.atoms``), so physics-based local stress clusters alongside
    the geometric/restraint problem types. One site per flagged atom -- including
    hydrogens; the spatial clustering turns a knot of flagged atoms (a distorted
    group + its H's + any clashing neighbour) into a single zone, while an isolated
    flag falls into the DBSCAN noise remainder.
    '''
    def __init__(self, atom, z):
        self._atom = atom
        self.z = z

    @property
    def atom(self):
        return self._atom

    @property
    def atoms(self):
        from chimerax.atomic import Atoms
        return Atoms([self._atom])

    @property
    def center(self):
        return self._atom.coord


class UnassessedResidueSite:
    '''
    A residue in the running simulation with NO local-strain reference -- i.e.
    none of its atoms has a baseline entry (a novel ligand, a modified residue,
    etc.). Surfaced as its own Problem-Zones category so that the strain tool's
    *silence* on such a residue is not mistaken for "no problem": we simply have
    no basis to judge it. Standard residues, including chain termini that carry a
    few unbaselined atoms (OXT, terminal H), are NOT flagged -- they still have
    assessed atoms.
    '''
    def __init__(self, residue):
        self._residue = residue

    @property
    def residue(self):
        return self._residue

    @property
    def atoms(self):
        return self._residue.atoms

    @property
    def center(self):
        return self._residue.atoms.coords.mean(axis=0)


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
            'Chiral outliers':                  self.get_chiral_problems,
            'Local strain':                     self.get_local_strain_problems,
            'Unassessed (strain)':              self.get_unassessed_strain_problems,
        }

    from chimerax.isolde.molobject import (
        ProperDihedralRestraint, DistanceRestraint,
        AdaptiveDihedralRestraint, AdaptiveDistanceRestraint,
        Rotamer, Rama, ChiralCenter
    )
    from chimerax.isolde.validation.clashes import Clash

    _registered_types = {
        'Standard torsion restraints': ProperDihedralRestraint,
        'Standard distance restraints': DistanceRestraint,
        'Adaptive torsion restraints': AdaptiveDihedralRestraint,
        'Reference distance restraints': AdaptiveDistanceRestraint,
        'Protein sidechains': Rotamer,
        'Protein backbone': Rama,
        'Clashes': Clash,
        'Chiral outliers': ChiralCenter,
        'Local strain': StrainSite,
        'Unassessed (strain)': UnassessedResidueSite,
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
    
    @staticmethod
    def _footprint_coords(site):
        '''
        The atomic footprint of a problem, used as its point set for spatial
        clustering. Every registered problem object exposes ``.atoms`` (the
        peptide unit for a Ramachandran outlier, the chi-dihedral sidechain for a
        rotamer, the pair for a clash, the flagged atom(s) for a strain hit),
        which gives DBSCAN a density that reflects the problem's real spatial
        extent rather than a single centroid. Falls back to ``.center`` for any
        exotic site lacking atoms.

        Note ``.atoms`` is not uniform across types: validation objects return an
        :class:`Atoms` collection, but the restraint objects return a plain tuple
        of :class:`Atom`. Coerce to :class:`Atoms` (as ``cluster_atoms`` does) so
        ``.coords`` is always available.
        '''
        import numpy
        from chimerax.atomic import Atoms
        atoms = getattr(site, 'atoms', None)
        if atoms is not None:
            if not isinstance(atoms, Atoms):
                atoms = Atoms(atoms)
            if len(atoms):
                return numpy.asarray(atoms.coords, dtype=float).reshape(-1, 3)
        center = getattr(site, 'center', None)
        if center is None:
            return None
        return numpy.asarray(center, dtype=float).reshape(1, 3)

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
                self.session.logger.warning(f'{validation_type} is not a type registered for problem aggregation!')
            else:
                sites.extend(vm(structure, outliers_only=validation_outliers_only))
        if not sites:
            return [], []
        # Densified clustering: each problem contributes its full atomic footprint
        # (see _footprint_coords) rather than a single centroid, so heterogeneous
        # problem types cluster on a common, spatially-honest density footing and a
        # contact-scale cutoff yields zones that hug the actual trouble instead of
        # ballooning. Points carry provenance back to their source problem.
        point_coords = []
        point_site = []
        for si, site in enumerate(sites):
            fp = self._footprint_coords(site)
            if fp is None or len(fp) == 0:
                continue
            point_coords.append(fp)
            point_site.extend([si] * len(fp))
        if not point_coords:
            return [], list(sites)
        coords = numpy.concatenate(point_coords, axis=0)
        point_site = numpy.array(point_site)
        from ..clustering import dbscan
        clusters, noise = dbscan(coords, cutoff, min_points)
        # Map each spatial cluster of points back to the unique source problems it
        # touches (a problem joins a zone if any of its footprint atoms fall in it).
        clustered_issues = []
        clustered_site_ids = set()
        for c in clusters:
            site_ids = numpy.unique(point_site[numpy.asarray(c, dtype=int)])
            clustered_site_ids.update(int(i) for i in site_ids)
            clustered_issues.append([sites[int(i)] for i in site_ids])
        remainder = [sites[i] for i in range(len(sites)) if i not in clustered_site_ids]
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

    @staticmethod
    def get_local_strain_problems(structure, outliers_only=True, z_threshold=None):
        '''
        Physics-based local stress problems, available only while `structure` is
        the model of a running ISOLDE simulation (the signal needs the live OpenMM
        System + current, near-minimised coordinates). Reconstructs the per-atom
        virial, subtracts the force-field frustration baseline, and flags atoms
        whose variance-normalised excess stress (von Mises z-score) exceeds a
        threshold. `outliers_only` picks the threshold per the Problem Zones
        convention (checked -> STRAIN_Z_OUTLIER ~3, unchecked -> STRAIN_Z_NONFAVORED
        ~2); pass `z_threshold` to override. Returns one StrainSite per flagged
        atom (or [] when no matching simulation is running).

        Coordinates are read from `all_atoms.coords` -- the snapshot ISOLDE
        already pushes to the GUI each update -- rather than from the live
        Context, so this is safe to call without pausing the (thread-stepped)
        simulation.
        '''
        session = structure.session
        isolde = getattr(session, 'isolde', None)
        if isolde is None or not isolde.simulation_running:
            return []
        sm = isolde.sim_manager
        sc = sm.sim_construct
        if sc.model is not structure:
            return []
        sh = sm.sim_handler
        system = getattr(sh, '_system', None)
        atoms = sc.all_atoms
        if system is None or len(atoms) != system.getNumParticles():
            return []
        import numpy
        from chimerax.isolde.validation.local_virial import LocalVirialCalculator
        from chimerax.isolde.validation.reference_stress import ReferenceStressLibrary
        pos_nm = atoms.coords / 10.0  # Angstrom -> nm, in System-particle order
        virial = LocalVirialCalculator(system).compute(pos_nm)['virial']
        try:
            lib = ReferenceStressLibrary()
        except FileNotFoundError:
            session.logger.warning('Local-strain problems: no reference stress '
                'library found; skipping. (Rebuild ISOLDE or run the builder.)')
            return []
        z = lib.compute_excess(virial, atoms)['von_mises_z']
        if z_threshold is None:
            z_threshold = STRAIN_Z_OUTLIER if outliers_only else STRAIN_Z_NONFAVORED
        flagged = numpy.where(numpy.nan_to_num(z, nan=-numpy.inf) > z_threshold)[0]
        return [StrainSite(atoms[int(i)], float(z[int(i)])) for i in flagged]

    @staticmethod
    def get_unassessed_strain_problems(structure, outliers_only=True):
        '''
        Residues in the running simulation for which the local-strain analysis has
        no reference at all (a residue is "unassessed" only if *none* of its atoms
        has a baseline entry -- novel ligands, modified residues). Emitted as its
        own category so the strain tool's silence on them is visible rather than
        read as "clean". Needs only the reference dictionary (no virial / coords),
        but is gated on an active simulation so it tracks exactly the atoms the
        Local strain stream would otherwise assess.
        '''
        session = structure.session
        isolde = getattr(session, 'isolde', None)
        if isolde is None or not isolde.simulation_running:
            return []
        sm = isolde.sim_manager
        sc = sm.sim_construct
        if sc.model is not structure:
            return []
        from chimerax.isolde.validation.reference_stress import ReferenceStressLibrary
        try:
            lib = ReferenceStressLibrary()
        except FileNotFoundError:
            return []
        sites = []
        for r in sc.all_atoms.unique_residues:
            if not any(lib.has_baseline(r.name, a.name) for a in r.atoms):
                sites.append(UnassessedResidueSite(r))
        return sites

    @staticmethod
    def get_chiral_problems(structure, outliers_only=True):
        # Chiral centres whose handedness is wrong or badly strained, judged from
        # the centre geometry (works on a static model, no simulation required).
        # Non-outlier centres are never "problems", so outliers_only is moot here.
        from chimerax.isolde.atomic.chirality import chiral_outliers
        chirals, oriented, severity = chiral_outliers(structure.session, structure.atoms)
        return chirals
