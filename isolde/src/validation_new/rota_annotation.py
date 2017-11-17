'''
Drawing code for live annotation of rotamers
'''
import numpy
from math import log

from chimerax.core.models import Drawing, Model
from chimerax.core.atomic import Bonds
from chimerax.core.geometry import translation


from ..color import standard_three_color_scale
from ..geometry import exclamation_mark, bond_cylinder_placements

class Rotamer_Annotations(Model):
    ''' Model holding annotations for current allowed or outlier rotamers. '''
    
    pickable = False
    
    def __init__(self, session, allowed_cutoff, outlier_cutoff):
        Model.__init__(self, 'Rotamer Annotations', session)
        
        # CA-CB bonds for all iffy (allowed or outlier) residues
        self._bonds = Bonds()
        self._scores = None
        self._current_colors = None
        
        from chimerax.core.atomic import get_triggers
        t = get_triggers(session)
        self._structure_change_handler = t.add_handler('changes', self.update_graphics)
        
        self._color_map = standard_three_color_scale('PiYG', log(outlier_cutoff), 0, log(allowed_cutoff))
        
        d = self._drawing = self._cb_annotation()
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
        if colors is None:
            # Update the colors on the next redraw
            self._update_needed = True
        else:
            self._current_colors = colors
    
    def update_graphics(self, *_):
        if self._scores is None or not len(self._scores):
            self.display = False
            return
        if not self.visible:
            return
        if self._update_needed:
            colors = self._current_colors = self._color_map.get_colors(numpy.log(self._scores))
        else:
            colors = self._current_colors
        bonds = self._bonds
        d = self._drawing
        d.positions = bond_cylinder_placements(bonds)
        d.colors = colors
        self._update_needed = False
        
        
        
    
    
    
    def _cb_annotation(self):
        from chimerax.core.surface.shapes import cylinder_geometry
        v, n, t = cylinder_geometry(radius=0.1, height=2, nc=8, caps=True)
        d = Drawing('rotamer cb indicator')
        d.vertices, d.normals, d.triangles = v, n, t
        return d
    
        
        
        
        
    
