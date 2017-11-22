'''
Drawing code for live annotation of rotamers
'''
import numpy
from math import log
from time import time

from chimerax.core.models import Drawing, Model
from chimerax.core.atomic import Bonds
from chimerax.core.geometry import translation, rotation, Places

from ..color import standard_three_color_scale
from ..geometry import exclamation_mark, spiral, bond_cylinder_placements
from ..geometry import scale_transforms

class Rotamer_Annotations(Model):
    ''' Model holding annotations for current allowed or outlier rotamers. '''
    
    pickable = False
    
    def __init__(self, session, atomic_model, allowed_cutoff, outlier_cutoff):
        Model.__init__(self, 'Rotamer Annotations', session)
        
        # CA-CB bonds for all iffy (allowed or outlier) residues
        self._bonds = Bonds()
        self._scores = None
        self._current_colors = None
        self._log_allowed_cutoff = log(allowed_cutoff)
        self._log_outlier_cutoff = log(outlier_cutoff)
        
        self._max_scale = 2
        
        from chimerax.core.atomic import get_triggers
        t = atomic_model.triggers
        self._structure_change_handler = t.add_handler('changes', self._update_graphics_if_needed)
        
        self._color_map = standard_three_color_scale('PiYG', log(outlier_cutoff), 0, log(allowed_cutoff))
        
        #d = self._drawing = self._cb_annotation()
        d = self._drawing = self._rota_indicator()
        self.add_drawing(d)
        self._update_needed = False
    
    @property
    def color_map(self):
        return self._color_map
    
    def delete(self):
        h = self._structure_change_handler
        if h is not None:
            from chimerax.core.atomic import get_triggers
            get_triggers(self.session).remove_handler(h)
            self._structure_change_handler = None
        Model.delete(self)
    
    def update_scores(self, iffy_bonds, iffy_scores, colors = None):
        if not len(iffy_bonds):
            self.display = False
            return
        self.display = True
        self._bonds = iffy_bonds
        self._scores = iffy_scores
        self._log_scores = numpy.log(iffy_scores)
        if colors is None:
            # Update the colors on the next redraw
            self._update_needed = True
        else:
            self._current_colors = colors
    
    def _update_graphics_if_needed(self, trigger_name, changes):
        changes = changes[1]
        if 'coord changed' in changes.atom_reasons():
            self.update_graphics()
    
    def update_graphics(self, *_, scale_by_scores = True):
        if self._scores is None or not len(self._scores):
            self.display = False
            return
        if not self.visible:
            return
        if self._update_needed:
            colors = self._current_colors = self._color_map.get_colors(self._log_scores)
        else:
            colors = self._current_colors
        bonds = self._bonds
        d = self._drawing
        if scale_by_scores:
            d.positions = self._scale_by_scores(bond_cylinder_placements(bonds).array())
        else:
            d.positions = bond_cylinder_placements(bonds)
        d.colors = colors
        self._update_needed = False
        
    def _scale_by_scores(self, transforms):
        log_scores = self._log_scores
        lac = self._log_allowed_cutoff
        loc = self._log_outlier_cutoff
        ms = self._max_scale
        scales = (log_scores-lac)/(loc-lac)*(ms-1)+1
        scales[scales > 2] = 2
        tf = scale_transforms(scales, transforms)
        return Places(place_array = tf)
        
        
        
        
    
    
    def _rota_indicator(self):
        v1, n1, t1 = exclamation_mark(radius=0.1, height=0.5, nc = 8)
        v2, n2, t2 = spiral(major_radius=0.3, minor_radius = 0.05, height=0.4,
                            turn_segments=6, circle_segments=3)
        translation((0,0,-0.15)).move(v2)
        v = numpy.concatenate((v1, v2))
        n = numpy.concatenate((n1, n2))
        t = numpy.concatenate((t1, t2+len(v1)))
        r = rotation((1,0,0),180)
        r.move(v)
        r.move(n)
        translation((0,0.5,0.25)).move(v)
        d = Drawing('rotamer indicator')
        d.vertices, d.normals, d.triangles = v, n, t
        return d
        
    
    def _exclamation_mark(self):
        v, n, t = exclamation_mark(radius=0.1, height=0.5, nc = 8)
        rotation((1,0,0),180).move(v)
        translation((0,1,0)).move(v)
        d = Drawing('rotamer cb indicator')
        d.vertices, d.normals, d.triangles = v, n, t
        return d


    
    def _cb_annotation(self):
        from chimerax.core.surface.shapes import cylinder_geometry
        v, n, t = cylinder_geometry(radius=0.1, height=2, nc=8, caps=True)
        d = Drawing('rotamer cb indicator')
        d.vertices, d.normals, d.triangles = v, n, t
        return d
    
        
        
        
        
    
