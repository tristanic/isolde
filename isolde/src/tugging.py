# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 12-Dec-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



from chimerax.core.models import Drawing, Model

from chimerax.mouse_modes import MouseMode
from .constants import defaults
_CARBON_MASS = defaults.CARBON_MASS

class TugAtomsMode(MouseMode):
    name = 'tug'
    #icon_file = 'tug.png'

    _modes = ('atom', 'residue', 'selection')

    def __init__(self, session, tug_mgr, atoms, spring_constant = None,
        mode = 'atom'):
        MouseMode.__init__(self, session)
        self.tug_mode = mode
        self._tugging = False
        self._last_xy = None
        self.name = 'ISOLDE_mouse_tug'
        self._focal_atom = None
        self._picked_atoms = None
        self._picked_tuggables = None
        self._pull_vector = None
        self._xyz0 = None
        self._xy = None
        # Atomic array to pick from
        self._atoms = atoms
        # Tuggable atoms manager
        self._tug_mgr = tug_mgr
        if spring_constant is None:
            from .constants import defaults
            spring_constant = defaults.MOUSE_TUG_SPRING_CONSTANT
        self.spring_constant = spring_constant
        self.structure = atoms.unique_structures[0]

        # Variables to be set by the caller
        self.last_tugged_atom = None
        self.already_tugging = False

    @property
    def tug_mode(self):
        return self._tug_mode

    @tug_mode.setter
    def tug_mode(self, mode):
        if mode not in self._modes:
            raise TypeError('Unsupported mode! Must be one of {}'.format(
                ','.join(self._modes)
            ))
        self._tug_mode = mode

    def __del__(self):
        self.cleanup()

    def cleanup(self):
        if self._picked_tuggables is not None:
            self._picked_tuggables.enableds = False

    @property
    def status(self):
        return self._tugging, self._picked_atom, self._xyz0
        print('Atom coords: ' + str(self._picked_atom.scene_coord) + ' xyz0: ' + str(self._xyz0))

    @property
    def tugging(self):
        return self._tugging

    @tugging.setter
    def tugging(self, flag):
        self._tugging = flag
        if not flag:
            self._focal_atom = None
            self._picked_atom = None
            if self._picked_tuggables is not None:
                self._picked_tuggables.enableds = False
                self._picked_tuggables = None

    def _pick_exclude(self, d):
        if not getattr(d, 'pickable', True):
            return True

        from chimerax.core.models import Model
        if isinstance(d, Model):
            if self.structure in d.all_models():
                return False
            return True
        if self.structure in d.drawing_lineage:
            return False
        return True

    def mouse_down(self, event):
        MouseMode.mouse_down(self, event)
        x,y = event.position()
        self._xy = (x,y)
        v = self.session.main_view
        if self.tug_mode == 'selection':
            pa = self._atoms[self._atoms.selecteds]
            self._reference_point = v.clip_plane_points(x,y)[0]
        else:
            # from . import picking
            # pick = picking.pick_closest_to_line(self.session, x, y, self._atoms, 0.5, hydrogens=True)
            pick = v.first_intercept(x, y, self._pick_exclude)
            pa = self._pick_atoms(pick)

        if pa is not None and len(pa):
            if self._start_tugging_atoms(pa):
                self._set_pull_direction(x, y)

    def _pick_atoms(self, pick):
        if pick is None:
            return None

        a = None
        r = None
        if hasattr(pick, 'atom'):
            a = pick.atom
            r = a.residue
        elif hasattr(pick, 'bond'):
            b = pick.bond
            coords = [a.coord for a in b.atoms]
            from chimerax.core.geometry import distance
            distances = [distance(c, pick.position) for c in coords]
            a = b.atoms[distances.index(min(distances))]
            r = a.residue
        elif hasattr(pick, 'residue'):
            r = pick.residue

        # Tug heavy atoms instead of hydrogens
        if a:
            if a.element.name=='H':
                h_mode = self._tug_mgr.allow_hydrogens
                if h_mode == 'no':
                    self.session.logger.warning('Tugging of hydrogens is not enabled. '
                        'Applying tug to the nearest bonded heavy atom.')
                    a = a.neighbors[0]
                elif h_mode == 'polar' and a.idatm_type == 'HC':
                    self.session.logger.warning('Tugging of non-polar hydrogens is not enabled. '
                        'Applying tug to the nearest bonded heavy atom.')
                    a = a.neighbors[0]
            self._focal_atom = a #a = self._focal_atom = pick

        tm = self.tug_mode
        if tm == "atom" and a:
            from chimerax.atomic import Atoms
            pa = Atoms([a])
        elif tm == "residue" and r:
            pa = r.atoms
            if a is not None:
                self._focal_atom = a
            else:
                self._focal_atom = r.principal_atom
        else:
            pa = None

        return pa

    def _start_tugging_atoms(self, atoms):
        tugs = self._picked_tuggables = self._tug_mgr.get_tuggables(atoms)
        n = len(tugs)
        if n == 0:
            self.tugging = False
            return False

        self._picked_atoms = pa = tugs.atoms
        tugs.targets = pa.coords

        # Scale the tugging force by atom masses and number of atoms
                # Scale spring constants down with atom count - we want
        # to be able to tug groups more strongly than single atoms, but not
        # *too* strongly.
        import numpy
        tugs.spring_constants = ((
            self.spring_constant * pa.elements.masses.astype(numpy.double)
            /_CARBON_MASS)/ n**(0.7)).reshape((n,1))
        tugs.enableds = True

        self.tugging = True
        return True

    def mouse_drag(self, event):
        if not self.tugging:
            return
        self._xy = x,y = event.position()
        ref_point = self._pull_reference_point()
        pull_vector = self._offset_vector(x, y, ref_point)
        tugs = self._picked_tuggables
        tugs.targets = self._picked_atoms.coords + pull_vector

    def mouse_up(self, event):
        MouseMode.mouse_up(self, event)
        self.tugging = False

    def _offset_vector(self, x, y, starting_coord):
        v = self.session.main_view
        x0,x1 = v.clip_plane_points(x, y)
        axyz = starting_coord
        # Project atom onto view ray to get displacement.
        dir = x1 - x0
        da = axyz - x0
        from chimerax.core.geometry import inner_product
        offset = self.structure.scene_position.inverse(is_orthonormal=True).transform_vectors(
            da - (inner_product(da, dir)/inner_product(dir,dir)) * dir
        )
        return -offset

    def _pull_reference_point(self):
        if self.tug_mode == 'selection':
            new_atom_center = self._picked_atoms.coords.mean(axis=0)
            lc = self._last_picked_atom_center
            if lc is not None:
                self._reference_point += (new_atom_center - lc)
            self._last_picked_atom_center = new_atom_center
            ref_point = self._reference_point
        else:
            ref_point = self._focal_atom.scene_coord
        return ref_point

    def _set_pull_direction(self, x, y):
        coords = self._picked_atoms.coords
        if self.tug_mode == 'selection':
            ref_point = self._reference_point
            self._last_picked_atom_center = coords.mean(axis=0)
        else:
            ref_point = self._focal_atom.scene_coord
        pull_vector = self._offset_vector(x, y, ref_point)
        tugs = self._picked_tuggables
        tugs.targets = coords + pull_vector

    def vr_press(self, event):
        # Virtual reality hand controller button press.
        if self.tug_mode == 'selection':
            pa = self._atoms[self._atoms.selecteds]
            self._reference_point = event.tip_position
            self._last_picked_atom_center = None
        else:
            view = self.session.main_view
            pick = event.picked_object(view)
            pa = self._pick_atoms(pick)

        if pa is not None and len(pa) > 0:
            self._start_tugging_atoms(pa)

    def vr_motion(self, event):
        # Virtual reality hand controller motion.
        if not self.tugging:
            return
        ref_point = self._pull_reference_point()
        pull_vector = event.tip_position - ref_point
        tugs = self._picked_tuggables
        tugs.targets = self._picked_atoms.coords + pull_vector
        # TODO: Apply torgue if there are multiple atoms.

    def vr_release(self, release):
        # Virtual reality hand controller button release.
        self.tugging = False



class HapticTugger():

    def __init__(self, session, index):
        self.session = session
        self._hh = session.HapticHandler
        self.index = index

        self.tugging = False

        self.tug_atom = None

    def start_tugging(self, atom):
        self.tug_atom = atom
        self.tugging = True
        self.update_target()
        #~ self._hh.setTargetPosition(self.index, *atom.coord, scene_coords = True)
        self._hh.startTugging(self.index)

    # Set the target position of the device to the current atom position, and
    # return the device position as the target for the atom
    def update_target(self):
        atom_xyz = self.tug_atom.coord
        pointer_xyz = self._hh.getPosition(self.index, scene_coords = True)
        self._hh.setTargetPosition(self.index, *atom_xyz, scene_coords = True)
        return pointer_xyz

    def stop_tugging(self):
        self._hh.stopTugging(self.index)
        self.tugging = False
        self.tug_atom = None

    def __del__(self):
        self.cleanup()

    def cleanup(self):
        if self.tugging:
            self.stop_tugging()
