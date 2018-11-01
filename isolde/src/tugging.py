# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



from chimerax.core.models import Drawing, Model

from chimerax.core.ui import MouseMode
from .constants import defaults
_CARBON_MASS = defaults.CARBON_MASS

class TugAtomsMode(MouseMode):
    name = 'tug'
    #icon_file = 'tug.png'

    _modes = ('atom', 'residue')

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

    def mouse_down(self, event):
        import numpy
        from chimerax.atomic import Atoms
        MouseMode.mouse_down(self, event)
        x,y = event.position()
        self._xy = (x,y)
        view = self.session.main_view
        from . import picking
        pick = picking.pick_closest_to_line(self.session, x, y, self._atoms, 0.5)
        if pick is not None:
            a = self._focal_atom = pick
            if self.tug_mode == "atom":
                pa = Atoms([a])
            else:
                pa = a.residue.atoms
            tugs = self._picked_tuggables = self._tug_mgr.get_tuggables(pa)
            pa = self._picked_atoms = tugs.atoms
            pull_vector = self._pull_direction(a, x, y)
            tugs.targets = tugs.atoms.coords + pull_vector
            # Scale the tugging force by atom masses and number of atoms
            n = len(tugs)
            tugs.spring_constants = ((
                self.spring_constant * pa.elements.masses.astype(numpy.double)
                /_CARBON_MASS)/ n).reshape((n,1))
            tugs.enableds = True

            self.tugging = True

    def mouse_drag(self, event):
        if not self.tugging:
            return
        self._xy = x,y = event.position()
        pull_vector = self._pull_direction(self._focal_atom, x, y)
        tugs = self._picked_tuggables
        tugs.targets = self._picked_atoms.coords + pull_vector

    def mouse_up(self, event):
        MouseMode.mouse_up(self, event)
        self.tugging = False

    def _pull_direction(self, atom, x, y):
        v = self.session.view
        x0,x1 = v.clip_plane_points(x, y)
        axyz = atom.scene_coord
        # Project atom onto view ray to get displacement.
        dir = x1 - x0
        da = axyz - x0
        from chimerax.core.geometry import inner_product
        offset = atom.structure.scene_position.inverse(is_orthonormal=True).transform_vectors(
            da - (inner_product(da, dir)/inner_product(dir,dir)) * dir
        )
        return -offset




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
