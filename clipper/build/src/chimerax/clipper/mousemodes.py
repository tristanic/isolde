import numpy
from chimerax.core.ui import mousemodes

class ZoomMouseMode(mousemodes.ZoomMouseMode):
    def __init__(self, session):
        mousemodes.ZoomMouseMode.__init__(self, session)
        self._far_clip = None
    
    @property
    def far_clip(self):
        if self._far_clip is None:
            v = self.session.view
            cp = v.clip_planes
            fc = cp.find_plane('far')
            if fc is not None:
                self._far_clip = fc
            else:
                from chimerax.core.graphics.view import ClipPlane
                c = v.camera
                cofr = v.center_of_rotation
                cpos = c.position.origin()
                clip_point = cpos + 2* (cofr-cpos)
                self._far_clip = ClipPlane('far', c.view_direction(), 
                                    clip_point, camera_normal = (0,0,1))
                v.clip_planes.add_plane(self._far_clip)
        return self._far_clip
                
    
    def zoom(self, delta_z):
        v = self.view
        c = v.camera
        if c.name() == 'orthographic':
            import numpy
            c.field_width = max(c.field_width - delta_z, self.pixel_size())
            # TODO: Make camera field_width a property so it knows to redraw.
            from chimerax.core.geometry import place
            cofr = v.center_of_rotation
            camera_to_cofr = cofr - c.position.origin()
            vd = c.view_direction()
            current_forward_distance = numpy.dot(camera_to_cofr, vd)
            new_view_distance = c.field_width * 0.9
            shift = current_forward_distance - new_view_distance
            new_origin = c.position.origin() + shift*vd
            new_pos = place.Place(axes = c.position.axes(), origin = new_origin)
            c.position = new_pos
            self.far_clip.plane_point = new_origin + (cofr-new_origin)*2
            c.redraw_needed = True
        else:
            shift = c.position.apply_without_translation((0, 0, delta_z))
            v.translate(shift)
            new_origin = c.position.origin()
            self.far_clip.plane_point = new_origin + (cofr-new_origin)*2
