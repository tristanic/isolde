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
            new_view_distance = c.field_width * 0.7
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


class SelectVolumeToContour(mousemodes.MouseMode):
    '''
    Designed to work together with ContourSelectedVolume.
    Each step of the mouse wheel increments through the currently 
    loaded Volume objects, temporarily selecting the current pick to
    highlight it in the display. Stores a reference to the last selected
    volume, accessible via picked_volume. If the last selected volume
    has never been set or has been deleted, picked_volume defaults to
    the first Volume object in session.models.list().
    '''
    def __init__(self, session):
        mousemodes.MouseMode.__init__(self, session)
        self._last_picked_index = 0
        self._picked_volume = None
        self._deselect_handler = None
        self._frames_until_deselect = 100
        self._deselect_counter = 0
    def wheel(self, event):
        d = int(event.wheel_value())
        vol_list = self._get_vol_list()
        for v in vol_list:
            v.selected = False
        n = len(vol_list)
        last = self._last_picked_index
        p = (last + d) % n
        sv = self._picked_volume = vol_list[p]
        sv.selected = True
        self._last_picked_index = p
        self.session.logger.status('Selected for contouring: {}'.format(sv.name))
        self._start_deselect_timer()
    def _get_vol_list(self):
        from chimerax.core.map import Volume
        return self.session.models.list(type=Volume)

    @property
    def picked_volume(self):
        try:
            vol_list = self._get_vol_list()
            if not len(vol_list):
                return None
            if self._picked_volume is None:
                self._picked_volume = vol_list[self._last_picked_index]
            elif self._picked_volume.deleted:
                self._picked_volume = vol_list[0]
                self._last_picked_index = 0
        except IndexError:
            self._picked_volume = vol_list[0]
            self._last_picked_index = 0
        return self._picked_volume
    
    def _start_deselect_timer(self):
        self._deselect_counter = 0
        if self._deselect_handler is None:
            self._deselect_handler = self.session.triggers.add_handler(\
                                'new frame', self._incr_deselect_counter)
    
    def _incr_deselect_counter(self, *_):
        self._deselect_counter += 1
        if self._deselect_counter >= self._frames_until_deselect:
            self.session.triggers.remove_handler(self._deselect_handler)
            self._deselect_handler = None
            self._picked_volume.selected = False

class ContourSelectedVolume(mousemodes.MouseMode):
    def __init__(self, session, selector, symmetrical=True):
        '''
        Modified volume contouring method which acts on a single volume
        at a time. By default, changes all contours towards/away from
        zero. If the volume has a surface_zone property set, it will be
        automatically masked back down after a short time delay.
        Args:
            selector:
                A SelectVolumeToContour object used to define the current
                target volume
            symmetrical:
                If True, scrolling up will adjust contours away from 
                zero (that is, negative contours will get more negative).
                If False, all contours will be shifted in the same
                direction.
        '''
        mousemodes.MouseMode.__init__(self, session)
        # SelectVolumeToContour object telling us which volume to work on
        self.selector = selector
        self.symmetrical = True
        self.target_volume = None
        self._remask_handler = None
        self._frames_until_remask = 10
        self._remask_counter = 0
    
    def wheel(self, event):
        d = event.wheel_value()
        v = self.selector.picked_volume
        if v is not None:
            self.target_volume = v
            if hasattr(v, 'overall_sigma'):
                sd = v.overall_sigma
            else:
                sd = v.mean_sd_rms()[1]
            step = d/30 * sd
            rep, levels = adjust_threshold_level(v, step, self.symmetrical)
            lsig = tuple(l/sd for l in levels)
            v.show()
            if rep != 'solid':
                lstr = ', '.join(format(l, '.3f') for l in levels)
                sstr = ', '.join(format(s, '.3f') for s in lsig)
                self.session.logger.status('Volume {} contour level(s): {} ({} sigma)'.format(v.name, lstr, sstr))
            self.session.ui.update_graphics_now()
            if hasattr(v, 'surface_zone'):
                if v.surface_zone.atoms is not None or v.surface_zone.coords is not None:
                    self._start_remask_countdown()
    
    def _start_remask_countdown(self):
        if self._remask_handler is None:
            self._remask_handler = self.session.triggers.add_handler('new frame', self._incr_remask_counter)
        self._remask_counter = 0
    
    def _incr_remask_counter(self, *_):
        self._remask_counter += 1
        if self._remask_counter >= self._frames_until_remask:
            from chimerax.core.surface.zone import surface_zone
            from chimerax.core
            v = self.target_volume
            if v.surface_zone.atoms is not None:
                coords = v.surface_zone.atoms.coords
                if v.surface_zone.coords is not None:
                    coords = numpy.concatenate([coords, v.surface_zone.coords])
            else:
                coords = v.surface_zone.coords
            
            surface_zone(v, coords, v.surface_zone.distance)
            self.session.triggers.remove_handler(self._remask_handler)
            self._remask_handler = None
            
    
def adjust_threshold_level(m, step, sym):
    if m.representation == 'solid':
        new_levels = [(l+step,b) for l,b in m.solid_levels]
        l,b = new_levels[-1]
        new_levels[-1] = (max(l,1.01*ms.maximum),b)
        m.set_parameters(solid_levels = new_levels)
    else:
        if sym and len(m.surface_levels) > 1:
            old_levels = m.surface_levels
            new_levels = []
            for l in old_levels:
                if l < 0:
                    newl = l-step
                    if newl > 0:
                        newl = -(abs(step))
                    new_levels.append(newl)
                else:
                    newl = l+step
                    if newl < 0:
                        newl = abs(step)
                    new_levels.append(newl)
        else:
            new_levels = tuple(l+step for l in m.surface_levels)
        m.set_parameters(surface_levels = new_levels)
    return(m.representation, new_levels)

        
        

    
    
