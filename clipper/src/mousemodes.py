# @Author: Tristan Croll
# @Date:   22-Mar-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 08-May-2019
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



import numpy
from chimerax.mouse_modes import MouseMode, ZoomMouseMode as ZoomMouseMode_Base

def initialize_clipper_mouse_modes(session):
    initialize_zoom_mouse_modes(session)
    initialize_map_contour_mouse_modes(session)

def initialize_zoom_mouse_modes(session):
    z = ZoomMouseMode(session)
    c = ClipPlaneAdjuster(session, z)
    session.ui.mouse_modes.bind_mouse_mode('right',['shift'], z)
    session.ui.mouse_modes.bind_mouse_mode('wheel',['shift'], c)

def initialize_map_contour_mouse_modes(session):
    #z = ZoomMouseMode(session)
    s = SelectVolumeToContour(session)
    v = ContourSelectedVolume(session, s, True)
    #session.ui.mouse_modes.bind_mouse_mode('right',[],z)
    session.ui.mouse_modes.bind_mouse_mode('wheel',['control'], s)
    session.ui.mouse_modes.bind_mouse_mode('wheel',[], v)

class ClipPlaneAdjuster(MouseMode):
    def __init__(self, session, zoom_mode):
        super().__init__(session)
        self._zoomer = zoom_mode

    def cleanup(self):
        pass

    def wheel(self, event):
        mult = 1-event.wheel_value()/30
        z = self._zoomer
        z.far_clip_multiplier *= mult
        z.near_clip_multiplier *= mult


class ZoomMouseMode(ZoomMouseMode_Base):
    def __init__(self, session):
        super().__init__(session)
        self._far_clip_multiplier = 0.5
        self._near_clip_multiplier = 0.5

    @property
    def far_clip_multiplier(self):
        '''
        Multiplier applied to the camera-cofr distance to decide the position
        of the rear clipping plane. Clamped to the range (0..1)
        '''
        return self._far_clip_multiplier

    @far_clip_multiplier.setter
    def far_clip_multiplier(self, val):
        val = max(min(val, 1), 0)
        self._far_clip_multiplier = val
        # Zoom with a delta-z of zero to force redraw
        self.zoom(0)

    @property
    def near_clip_multiplier(self):
        '''
        Multiplier applied to the camera-cofr distance to decide the position
        of the near clipping plane. Clamped to the range (0..1)
        '''
        return self._near_clip_multiplier

    @near_clip_multiplier.setter
    def near_clip_multiplier(self, val):
        val = max(min(val, 1), 0)
        self._near_clip_multiplier = val
        # Zoom with a delta-z of zero to force redraw
        self.zoom(0)

    @property
    def far_clip(self):
        v = self.session.view
        cp = v.clip_planes
        fc = cp.find_plane('far')
        if fc is not None:
            return fc
        else:
            from chimerax.core.graphics.view import ClipPlane
            c = v.camera
            cofr = v.center_of_rotation
            cpos = c.position.origin()
            clip_point = cpos + (1+self.far_clip_multiplier)* (cofr-cpos)
            fc = ClipPlane('far', c.view_direction(),
                                clip_point, camera_normal = (0,0,1))
            # Put the near clip at the camera position for depth cueing
            v.clip_planes.add_plane(fc)
            return fc

    @property
    def near_clip(self):
        v = self.session.view
        cp = v.clip_planes
        nc = cp.find_plane('near')
        if nc is not None:
            return nc
        else:
            from chimerax.core.graphics.view import ClipPlane
            c = v.camera
            cofr = v.center_of_rotation
            cpos = c.position.origin()
            clip_point = cpos + (1-self.near_clip_multiplier) * (cofr-cpos)
            nc = ClipPlane('near', c.view_direction(),
                                clip_point, camera_normal = (0,0,-1))
            # Put the near clip at the camera position for depth cueing
            v.clip_planes.add_plane(nc)
            return nc

    def cleanup(self):
        pass

    def zoom(self, delta_z, stereo_scaling = False):
        v = self.view
        c = v.camera
        cofr = v.center_of_rotation
        if stereo_scaling and c.name == 'stereo':
            v.stereo_scaling(delta_z)
        if c.name == 'orthographic':
            import numpy
            c.field_width = max(c.field_width - delta_z, self.pixel_size())
            # TODO: Make camera field_width a property so it knows to redraw.
            from chimerax.core.geometry import place
            camera_to_cofr = cofr - c.position.origin()
            vd = c.view_direction()
            current_forward_distance = numpy.dot(camera_to_cofr, vd)
            new_view_distance = c.field_width * 0.7
            shift = current_forward_distance - new_view_distance
            new_origin = c.position.origin() + shift*vd
            new_pos = place.Place(axes = c.position.axes(), origin = new_origin)
            c.position = new_pos
            distance = cofr-new_origin
            self.far_clip.plane_point = new_origin + distance*(1+self.far_clip_multiplier)
            self.near_clip.plane_point = new_origin + distance*(1-self.near_clip_multiplier)
            c.redraw_needed = True
        else:
            shift = c.position.transform_vectors((0, 0, delta_z))
            v.translate(shift)
            new_origin = c.position.origin()
            self.far_clip.plane_point = new_origin + (cofr-new_origin)*2

class SelectVolumeToContour(MouseMode):
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
        super().__init__(session)
        self._last_picked_index = 0
        self._picked_volume = None
        self._deselect_handler = None
        self._time_until_deselect = 5 # seconds
        self._deselect_start_time = 0
    def wheel(self, event):
        '''Select the next visible volume.'''
        d = int(event.wheel_value())
        vol_list = self._get_vol_list()
        for v in vol_list:
            v.selected = False
        n = len(vol_list)
        if n == 0:
            return
        last = self._last_picked_index
        p = (last + d) % n
        sv = self._picked_volume = vol_list[p]
        sv.selected = True
        self._last_picked_index = p
        self.session.logger.status('Selected for contouring: {}'.format(sv.name))
        self._start_deselect_timer()
    def _get_vol_list(self):
        '''Get a list of currently visible volumes.'''
        from chimerax.map import Volume
        vlist = self.session.models.list(type=Volume)
        for i in reversed(range(len(vlist))):
            if not vlist[i].visible:
                vlist.pop(i)
        return vlist


    @property
    def picked_volume(self):
        try:
            vol_list = self._get_vol_list()
            if not len(vol_list):
                return None
            pv = self._picked_volume
            if pv is None:
                pv = self._picked_volume = vol_list[self._last_picked_index]
            elif pv not in vol_list or not pv.visible:
                pv = self._picked_volume = vol_list[0]
                self._last_picked_index = 0
        except IndexError:
            pv = self._picked_volume = vol_list[0]
            self._last_picked_index = 0
        return self._picked_volume

    def _start_deselect_timer(self):
        from time import time
        self._deselect_start_time = time()
        if self._deselect_handler is None:
            self._deselect_handler = self.session.triggers.add_handler(\
                                'new frame', self._deselect_on_timeout)

    def _deselect_on_timeout(self, *_):
        from time import time
        if time()- self._deselect_start_time > self._time_until_deselect:
            self.session.triggers.remove_handler(self._deselect_handler)
            self._deselect_handler = None
            self._picked_volume.selected = False

class ContourSelectedVolume(MouseMode):
    def __init__(self, session, selector, symmetrical=True):
        '''
        Modified volume contouring method which acts on a single volume
        at a time. By default, changes all contours towards/away from
        zero. If the volume has a surface_zone property set, it will be
        automatically masked back down after contouring.
        Args:
            session:
                The ChimeraX session.
            selector:
                A SelectVolumeToContour object used to define the current
                target volume.
            symmetrical:
                If True, scrolling up will adjust contours away from
                zero (that is, negative contours will get more negative).
                If False, all contours will be shifted in the same
                direction.
        '''
        super().__init__(session)
        # SelectVolumeToContour object telling us which volume to work on
        self.selector = selector
        self.symmetrical = symmetrical
        self.target_volume = None


    def wheel(self, event):
        d = event.wheel_value()
        v = self.selector.picked_volume
        if v is not None:
            self.target_volume = v
            if hasattr(v, 'sigma'):
                sd = v.sigma
            else:
                sd = v.mean_sd_rms()[1]
            step = d/30 * sd
            rep, levels = adjust_threshold_level(v, step, self.symmetrical)
            lsig = tuple(l/sd for l in levels)
            if rep != 'solid':
                lstr = ', '.join(format(l, '.3f') for l in levels)
                sstr = ', '.join(format(s, '.3f') for s in lsig)
                self.session.logger.status('Volume {} contour level(s): {} ({} sigma)'.format(v.name, lstr, sstr))




def adjust_threshold_level(m, step, sym):
    if m.surfaces_in_style('solid'):
        new_levels = [(l+step,b) for l,b in m.solid_levels]
        l,b = new_levels[-1]
        new_levels[-1] = (max(l,1.01*ms.maximum),b)
        m.set_parameters(solid_levels = new_levels)
        return ('solid', new_levels)
    else:
        #if sym and len(m.surface_levels) > 1:
        if sym and len(m.surfaces) > 1:
            #old_levels = m.surface_levels
            old_levels = [s.level for s in m.surfaces]
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
            old_levels = [s.level for s in m.surfaces]
            #new_levels = tuple(l+step for l in m.surface_levels)
            new_levels = tuple(l+step for l in old_levels)
        m.set_parameters(surface_levels = new_levels)
    return (None, new_levels)
