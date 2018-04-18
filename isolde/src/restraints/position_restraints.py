# @Author: Tristan Croll
# @Date:   14-Feb-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   Tristan Croll
# @Last modified time: 18-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



import numpy

from chimerax.core.models import Model, Drawing
from chimerax.core.geometry import Places

from ..molobject import Position_Restraint_Mgr, Position_Restraint
from ..molarray import Position_Restraints


class Live_Position_Restraint_Mgr(Position_Restraint_Mgr):
    '''
    Manages position restraints and their visualisations for a single atomic
    structure.
    '''
    def __init__(self, model, c_pointer=None):
        super().__init__(model, c_pointer)
        self._visibles = None
        pd = self._pin_drawing = Drawing('Target pins')
        bd = self._bond_drawing = Drawing('Restraint bonds')
        self.add_drawing(pd)
        self.add_drawing(bd)
        pd.vertices, pd.normals, pd.triangles = _target_pin_geometry()
        bd.vertices, bd.normals, bd.triangles = _pseudobond_geometry()
        self.set_bond_color(_default_bond_color)
        self.set_pin_color(_default_pin_color)
        pd.display = False
        bd.display = False
        self.model.triggers.add_handler('changes', self._update_graphics_if_needed)

    def _update_visibles(self):
        self._visibles = self.visible_restraints

    @property
    def visibles(self):
        if self._visibles is None:
            self._update_visibles()
        return self._visibles

    def _update_graphics_if_needed(self, trigger_name, changes):
        update_needed = False
        atom_reasons = changes[1].atom_reasons()
        if 'display changed' in atom_reasons or 'hide changed' in atom_reasons:
            self._update_visibles()
            update_needed = True
        if 'coord changed' in atom_reasons:
            update_needed = True
        if update_needed:
            self._update_graphics()

    def _update_graphics(self):
        pd = self._pin_drawing
        bd = self._bond_drawing
        visibles = self.visibles
        n = len(visibles)
        if n==0:
            pd.display = False
            bd.display = False
            return
        pd.display = True
        bd.display = True
        self._update_pin_drawing(pd, visibles, n)
        self._update_bond_drawing(bd, visibles, n)

    def _update_pin_drawing(self, pd, visibles, n):
        xyzr = numpy.ones((n,4), numpy.double)
        xyzr[:, :3] = visibles.targets
        pd.positions = Places(shift_and_scale=xyzr)

    def _update_bond_drawing(self, bd, visibles, n):
        bd.positions = Places(opengl_array = visibles._bond_cylinder_transforms)

    def set_pin_color(self, color):
        self._pin_drawing.color = color

    def set_bond_color(self, color):
        self._bond_drawing.color = color







_default_bond_color = [200, 250, 160, 255]
_default_pin_color = [200,0,200,255]

def _pseudobond_geometry(segments=9):
    from chimerax.core import surface
    return surface.dashed_cylinder_geometry(segments, height=0.5)

def _target_pin_geometry(handle_radius=0.4, pin_radius=0.05, total_height=1.0):
    from ..geometry import pin_geometry
    return pin_geometry(handle_radius, pin_radius, total_height)
