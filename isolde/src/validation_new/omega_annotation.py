
import numpy

from chimerax.core.models import Model, Drawing

from .. import geometry

class Omega_Annotation(Model):
    def __init__(self, session, parent,
                    non_pro_cis_color = [255,32,32,255],
                    pro_cis_color = [32,255,32,255],
                    twisted_color = [255,255,32,255]):
        '''
        Handles drawing of the "filled-in cup" annotation for cis and
        twisted peptide bonds.
        '''
        self._non_pro_cis_color = numpy.array(non_pro_cis_color).astype('uint8')
        self._pro_cis_color = numpy.array(pro_cis_color).astype('uint8')
        self._twisted_color = numpy.array(twisted_color).astype('uint8')
        super().__init__('cis/twisted peptides', session)
        d = self._drawing = Drawing('omega annotations')
        self.add_drawing(d)
        parent.add([self])
    
    def update_coords(self, dihedrals, twisted_mask, cis_pro_mask):
        '''
        Update the drawing.
        Args:
            dihedrals:
                a Dihedrals object containing the set of problem dihedrals
            twisted_indices:
                indices for all vertices involved in a twisted peptide
                bond
            cis_pro_indices:
                indices for all vertices involved in a cis-proline 
                peptide bond
        
        All other vertices will be coloured with the non-pro cis colour.
        '''
        d = self._drawing
        geometry.dihedral_fill_and_color_planes(dihedrals, d, 
            twisted_mask, cis_pro_mask, self._non_pro_cis_color,
            self._twisted_color, self._pro_cis_color)
        
    
       
