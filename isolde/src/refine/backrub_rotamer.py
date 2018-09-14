
class Backrub:
    '''
    Automatically fit a protein sidechain into density using the
    Backrub algorithm (Davis et al. 2006, Structure 14: 265-274)
    '''
    def __init__(self, residue, density_map):
        self._density_map = density_map
        from chimerax.atomic import Residue
        if residue.polymer_type != Residue.PT_AMINO:
            raise TypeError('Backrub is only applicable to protein residues!')
        from chimerax.isolde import session_extensions as sx
        rota_mgr = sx.get_rotamer_mgr(residue.session)
        rota = self.rotamer = rota_mgr.get_rotamer(residue)
        if rota is None:
            raise TypeError('Not a rotameric residue!')
        res_atoms = self._main_res_atoms = residue.atoms
        natom = res_atoms[res_atoms.names=='N'][0]
        catom = res_atoms[res_atoms.names=='C'][0]
        neighbors = residue.neighbors
        prev_res = None
        next_res = None
        if len(neighbors) > 2:
            raise TypeError('Too many residues bonded to the target! Backrub will not work'
                ' with post-translational modifications.')
        for neighbor in neighbors:
            bond = residue.bonds_between(neighbor)[0]
            batoms = bond.atoms
            if natom in batoms:
                prev_res = neighbor
            elif catom in batoms:
                next_res = neighbor
        if prev_res is None or next_res is None:
            raise TypeError('Backrub requires N- and C-terminal flanking residues to be present!')

        prev_atoms = self._prev_res_atoms = prev_res.atoms
        next_atoms = self._next_res_atoms = next_res.atoms

        prev_ca = prev_atoms[prev_atoms.names=='CA'][0]
        next_ca = next_atoms[next_atoms.names=='CA'][0]
        this_ca = residue.atoms[residue.atoms.names=='CA'][0]
        import numpy
        self.primary_axis = prev_ca.coord - next_ca.coord
        self._primary_axis_center = (prev_ca.coord + next_ca.coord)/2
        self.prev_pep_axis = this_ca.coord - prev_ca.coord
        self._prev_pep_axis_center = (this_ca.coord + prev_ca.coord)/2
        self.next_pep_axis = next_ca.coord - this_ca.coord
        self._next_pep_axis_center = (next_ca.coord + this_ca.coord)/2

        self._all_movable_atoms = res_atoms.merge(
            prev_atoms[numpy.in1d(prev_atoms.names, ['C', '0'])]
        ).merge(
            next_atoms[numpy.in1d(next_atoms.names, ['N', 'H'])]
        )

        self._prev_pep = res_atoms[numpy.in1d(res_atoms.names, ['N', 'H'])].merge(
            prev_atoms[numpy.in1d(prev_atoms.names, ['C', 'O'])]
        )

        self._next_pep = res_atoms[numpy.in1d(res_atoms.names, ['C', 'O'])].merge(
            prev_atoms[numpy.in1d(prev_atoms.names, ['N', 'H'])]
        )

    def auto_fit(self):
        '''
        Attempt to automatically find the rotamer that best fits the density.
        '''
        # Store the original coordinates to allow reversion
        self._original_coords = self._all_movable_atoms.coords
        import numpy
        from scipy.optimize import minimize
        # First optimize the rotation about the principal axis to get the best
        # fit for CA and CB
        res_atoms = self._main_res_atoms
        moving_atoms = self._all_movable_atoms
        check_atoms = res_atoms[numpy.in1d(res_atoms.names, ['CA', 'CB'])]
        minimize(self.rotate_and_check_fit, 0,
            args=(self.primary_axis, self._primary_axis_center, moving_atoms, check_atoms),
            bounds=[[-10,10]])

        moving_atoms = self._prev_pep
        check_atoms = self._prev_pep
        minimize(self.rotate_and_check_fit, 0,
            args=(self.prev_pep_axis, self._prev_pep_axis_center, moving_atoms, check_atoms),
            bounds=[[-10,10]])

        moving_atoms = self._next_pep
        check_atoms = self._next_pep
        minimize(self.rotate_and_check_fit, 0,
            args=(self.next_pep_axis, self._next_pep_axis_center, moving_atoms, check_atoms),
            bounds=[[-10,10]])

        rotamer = self.rotamer
        num_targets = rotamer.num_targets
        targets = [rotamer.get_target(i) for i in range(num_targets)]
        num_chi = rotamer.num_chi_dihedrals
        current_chi = 0

        cutoff_frac = 0.9

        # from chimerax.atomic import Atoms
        # from math import degrees
        # for i in range(num_chi):
        #     chi = rotamer.chi_dihedrals[i]
        #     moving_atoms = rotamer.moving_atoms(i)
        #     check_atoms = []
        #     for b in chi.atoms[2].bonds:
        #         check_atoms.extend(b.atoms)
        #     check_atoms = Atoms(check_atoms)
        #     results = []
        #     coords = chi.atoms.coords
        #     axis = coords[2]-coords[1]
        #     from chimerax.core.geometry import matrix
        #     center = matrix.project_to_axis(coords[3], axis, coords[1])
        #     for t in targets:
        #         results.append(self.rotate_and_check_fit(
        #             [degrees(t['Angles'][i]-chi.angle)], axis, center,
        #             moving_atoms, check_atoms
        #         ))
        #     print("Results for chi {}: {}".format(i, results))
        #     min_val = min(results)
        #     cutoff_val = min_val*cutoff_frac
        #     new_targets = []
        #     for t, r in zip(targets, results):
        #         if r < cutoff_val:
        #             new_targets.append(t)
        #     targets = new_targets
        #     if len(targets) == 1: break

        # By now we should have narrowed it down to quite a short list
        # (hopefully just one). Now to do some more detailed minimization
        print("Number of remaining targets: {}".format(len(targets)))
        print("Remaining targets: {}".format(
            ','.join([t['Name'] for t in targets])
        ))

        final_results = []
        for target in targets:
            t_angles = target['Angles']
            esds = target['ESDs']
            bounds = [[t-e, t+e] for t,e in zip(t_angles, esds)]
            print('minimizing target {}'.format(target['Name']))
            print('Target angles: {}'.format(t_angles))
            print('Bounds: {}'.format(bounds))
            final_results.append(
                minimize(self.fine_tune_rotamer, t_angles,
                args=(rotamer,), bounds=bounds ))
        print(len(final_results))
        scores = numpy.array([r.fun for r in final_results])
        print("Scores: {}".format(scores))
        best = final_results[numpy.argmin(scores)].x
        self.fine_tune_rotamer(best, rotamer)





    def rotate_and_check_fit(self, angle, axis, center, moving_atoms, check_atoms):
        '''
        Returns the negative sum of density values at the centre of each of the atoms in
        check_atoms after rotating moving_atoms by angle around axis.
        '''
        angle = angle[0]
        from chimerax.core.geometry import Place, rotation
        coords = moving_atoms.coords
        tf = rotation(axis, angle, center)
        moving_atoms.coords = tf.moved(moving_atoms.coords)
        result = -sum(self._density_map.interpolated_values(check_atoms.coords))
        return result

    def fine_tune_rotamer(self, angles, rotamer):
        check_atoms = rotamer.residue.atoms
        import numpy
        from chimerax.core.geometry import Place, rotation, matrix
        #angles = numpy.array(angles)
        rot_angles = numpy.degrees(angles-rotamer.angles)
        for i in range(rotamer.num_chi_dihedrals):
            ma = rotamer.moving_atoms(i)
            chi = rotamer.chi_dihedrals[i]
            coords = chi.atoms.coords
            axis = coords[2]-coords[1]
            center = matrix.project_to_axis(coords[3], axis, coords[1])
            tf = rotation(axis, rot_angles[i], center)
            ma.coords = tf.moved(ma.coords)
        result = -sum(self._density_map.interpolated_values(check_atoms.coords))
        print("Result: {}".format(result))
        return result
        # return -sum(self._density_map.interpolated_values(check_atoms.coords))

def test_backrub(session, isolde):
    m = isolde.selected_model
    r = m.residues[70]
    from ..view import focus_on_selection
    focus_on_selection(session, r.atoms)
    from chimerax.clipper.symmetry import get_symmetry_handler
    sh = get_symmetry_handler(m)
    mdff_p = sh.xmapset['MDFF potential']
    b = Backrub(r, mdff_p)
    b.auto_fit()
