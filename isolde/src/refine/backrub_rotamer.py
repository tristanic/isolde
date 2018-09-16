
class Backrub:
    '''
    Automatically fit a protein sidechain into density using the
    Backrub algorithm (Davis et al. 2006, Structure 14: 265-274)
    '''
    def __init__(self, residue, density_map):
        dm = self._density_map = density_map
        # Make sure the map is covering the selected residue
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
            prev_atoms[numpy.in1d(prev_atoms.names, ['C', 'O'])]
        ).merge(
            next_atoms[numpy.in1d(next_atoms.names, ['N', 'H'])]
        )

        self._prev_pep = res_atoms[numpy.in1d(res_atoms.names, ['N', 'H'])].merge(
            prev_atoms[numpy.in1d(prev_atoms.names, ['C', 'O'])]
        )

        self._next_pep = res_atoms[numpy.in1d(res_atoms.names, ['C', 'O'])].merge(
            next_atoms[numpy.in1d(next_atoms.names, ['N', 'H'])]
        )

    def auto_fit(self):
        '''
        Attempt to automatically find the rotamer that best fits the density.
        '''
        dm = self._density_map
        if not dm.display:
            # Make sure the map is updated to cover the current region
            dm._box_changed_cb('auto_fit', (dm.box_params, True))
        map_range = dm.stats.max - dm.stats.min
        # Store the original coordinates to allow reversion
        self._original_coords = self._all_movable_atoms.coords
        import numpy
        from scipy.optimize import minimize, minimize_scalar
        # First optimize the rotation about the principal axis to get the best
        # fit for CA and CB
        res_atoms = self._main_res_atoms
        moving_atoms = self._all_movable_atoms
        check_atoms = res_atoms[numpy.in1d(res_atoms.names, ['CA', 'CB'])]
        original_coords = moving_atoms.coords
        #self.rotate_and_check_fit([0], self.primary_axis, self._primary_axis_center, moving_atoms, original_coords, check_atoms)
        minimize_scalar(self.rotate_and_check_fit, bounds=[-20,20],
            args=(self.primary_axis, self._primary_axis_center, moving_atoms, original_coords, check_atoms),
            )

        moving_atoms = self._prev_pep
        check_atoms = self._prev_pep
        original_coords = moving_atoms.coords
        minimize_scalar(self.rotate_and_check_fit, bounds=[-20,20],
            args=(self.prev_pep_axis, self._prev_pep_axis_center, moving_atoms, original_coords, check_atoms),
            )

        moving_atoms = self._next_pep
        check_atoms = self._next_pep
        original_coords = moving_atoms.coords
        minimize_scalar(self.rotate_and_check_fit, bounds=[-20,20],
            args=(self.next_pep_axis, self._next_pep_axis_center, moving_atoms, original_coords, check_atoms),
            )

        rotamer = self.rotamer
        num_targets = rotamer.num_targets
        targets = [rotamer.get_target(i) for i in range(num_targets)]
        num_chi = rotamer.num_chi_dihedrals
        current_chi = 0

        cutoff_frac = 0.9

        from chimerax.atomic import Atoms
        from math import degrees
        for i in range(num_chi):
            chi = rotamer.chi_dihedrals[i]
            moving_atoms = rotamer.moving_atoms(i)
            check_atoms = []
            for b in chi.atoms[2].bonds:
                check_atoms.extend(b.atoms)
            check_atoms = Atoms(check_atoms)
            results = []
            for t in targets:
                # rotate all the previously checked chis to this rotamer
                for j in range(i):
                    c = rotamer.chi_dihedrals[j]
                    ma = rotamer.moving_atoms(j)
                    self.rotate_chi_to(c, ma, t['Angles'][j])


                results.append(self.rotate_chi_and_check_fit(
                    chi, moving_atoms, t['Angles'][i], check_atoms
                ))
            results = numpy.array(results)
            # print("Results for chi {}: {}".format(i, results))
            # If there is a wide range of values, select out the top few.
            # Otherwise, keep them all
            if max(results) - min(results) > map_range/10:
                min_val = min(results)
                # Normalise to the range 0..1, and discard everything more than 20%
                # away from the minimum
                results -= min_val
                results /= max(results)
                cutoff_val = 0.2
                new_targets = []
                for t, r in zip(targets, results):
                    if r < cutoff_val:
                        new_targets.append(t)
                targets = new_targets
            if len(targets) == 1: break

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
            # print('minimizing target {}'.format(target['Name']))
            # print('Target angles: {}'.format(t_angles))
            # print('Bounds: {}'.format(bounds))
            final_results.append(
                minimize(self.fine_tune_rotamer, t_angles,
                args=(rotamer,), bounds=bounds ))
        #print(len(final_results))
        scores = numpy.array([r.fun for r in final_results])
        # print("Scores: {}".format(scores))
        best = final_results[numpy.argmin(scores)].x
        self.fine_tune_rotamer(best, rotamer)




    def rotate_chi_to(self, chi_dihedral, moving_atoms, target_angle):
        ma = moving_atoms
        from math import degrees
        coords = chi_dihedral.atoms.coords
        axis = coords[2]-coords[1]
        from chimerax.core.geometry import matrix
        center = matrix.project_to_axis(coords[3], axis, coords[1])
        from chimerax.core.geometry import rotation
        tf = rotation(axis, degrees(target_angle-chi_dihedral.angle), center)
        ma.coords = tf.moved(ma.coords)

    def rotate_chi_and_check_fit(self, chi_dihedral, moving_atoms, target_angle, check_atoms = None):
        '''
        Returns the negative sum of density values at the centre of each of the atoms in
        check_atoms after rotating moving_atoms by angle around axis.
        '''
        if check_atoms is None:
            check_atoms = moving_atoms
        weights = check_atoms.elements.numbers
        self.rotate_chi_to(chi_dihedral, moving_atoms, target_angle)
        result = -sum(self._density_map.interpolated_values(check_atoms.coords)*weights)/sum(weights)
        return result

    def rotate_and_check_fit(self, angle, axis, center, moving_atoms, original_coords, check_atoms):
        '''
        Returns the negative sum of density values at the centre of each of the atoms in
        check_atoms after rotating moving_atoms by angle around axis.
        '''
        from chimerax.core.geometry import Place, rotation
        coords = moving_atoms.coords
        weights = check_atoms.elements.numbers
        tf = rotation(axis, angle, center)
        moving_atoms.coords = tf.moved(original_coords)
        dvals, outside = self._density_map.interpolated_values(check_atoms.coords, out_of_bounds_list=True)
        if len(outside):
            raise RuntimeError('At least one atom is currently projecting past'
                ' the edge of the displayed map box. Re-center the map on the'
                ' residue before trying again')
        result = -sum(dvals*weights)/sum(weights)
        print("Result for angle {}: {}".format(angle, result))
        return result
        #return -sum(dvals)


    def fine_tune_rotamer(self, angles, rotamer):
        check_atoms = rotamer.residue.atoms
        nchi = rotamer.num_chi_dihedrals
        import numpy
        weights = check_atoms.elements.numbers
        from chimerax.core.geometry import Place, rotation, matrix
        #angles = numpy.array(angles)
        rot_angles = numpy.degrees(angles-rotamer.angles)
        for i in range(nchi):
            ma = rotamer.moving_atoms(i)
            weights[check_atoms.indices(ma)] *= 2**(nchi-i)
            chi = rotamer.chi_dihedrals[i]
            coords = chi.atoms.coords
            axis = coords[2]-coords[1]
            center = matrix.project_to_axis(coords[3], axis, coords[1])
            tf = rotation(axis, rot_angles[i], center)
            ma.coords = tf.moved(ma.coords)
        dvals, outside = self._density_map.interpolated_values(check_atoms.coords,
            out_of_bounds_list = True)
        if len(outside):
            raise RuntimeError('At least one atom is currently projecting past'
                ' the edge of the displayed map box. Re-center the map on the'
                ' residue before trying again')
        result = -sum(dvals*weights)/sum(weights)
        return result
        # return -sum(self._density_map.interpolated_values(check_atoms.coords))

def _true():
    return True

def test_backrub(session, isolde, selected=False, focus=True):
    m = isolde.selected_model
    if selected:
        from chimerax.atomic import selected_atoms
        r = selected_atoms(session).unique_residues[0]
    else:
        r = m.residues[36]
    if focus:
        from ..view import focus_on_selection
        focus_on_selection(session, r.atoms)
    from chimerax.clipper.symmetry import get_symmetry_handler
    sh = get_symmetry_handler(m)
    mdff_p = sh.xmapset['MDFF potential']
    b = Backrub(r, mdff_p)
    from ..delayed_reaction import delayed_reaction
    delayed_reaction(session.triggers, 'new frame', _true, [], _true, b.auto_fit, [])
    #b.auto_fit()
