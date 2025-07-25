# @Author: Tristan Croll <tic20>
# @Date:   26-Apr-2020
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 24-Nov-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def get_stepper_mgr(session):
    if hasattr(session, '_isolde_steppers'):
        return session._isolde_steppers
    return ResidueStepperMgr(session)

def get_stepper(structure, session_restore=False):
    '''
    Get the :class:`ResidueStepper` controlling ISOLDE's navigation around the
    given structure, creating it if it doesn't yet exist.
    '''    
    return get_stepper_mgr(structure.session).get_stepper(structure)

from chimerax.core.state import StateManager, State

class ResidueStepperMgr(StateManager):
    def __init__(self, session):
        self._steppers = {}
        session._isolde_steppers = self
        if not len(session.state_managers(ResidueStepperMgr)):
            self.init_state_manager(session, base_tag='Isolde Residue Stepper Manager')
    
    def register_stepper(self, stepper):
        m = stepper.structure
        if stepper.structure in self._steppers.keys():
            raise KeyError(f'Model #{m.id_string} already has a stepper registered!')
        self._steppers[m] = stepper
    
    def remove_stepper(self, stepper):
        self._steppers.pop(stepper.structure, None)
    
    def get_stepper(self, structure):
        if structure in self._steppers.keys():
            return self._steppers[structure]
        return ResidueStepper(structure)
    
    def take_snapshot(self, session, flags):
        from . import ISOLDE_STATE_VERSION
        data = {
            'version': ISOLDE_STATE_VERSION,
            'steppers': self._steppers
        }
        return data
    
    def reset_state(self, session):
        self._steppers = {}

    @staticmethod
    def restore_snapshot(session, data):
        mgr = ResidueStepperMgr(session)
        mgr._steppers = data['steppers']
        return mgr




    


class ResidueStepper(State):

    _num_created = 0

    DEFAULT_INTERPOLATE_FRAMES=15
    DEFAULT_MAX_INTERPOLATE_DISTANCE=10
    DEFAULT_VIEW_DISTANCE=10
    DEFAULT_DIRECTION="next"
    '''
    Provide methods to step forward/backward through the polymeric residues in
    a structure.
    '''
    def __init__(self, structure, view_distance=12, session_restore=False):
        self.session = structure.session
        self.structure = structure
        self._current_residue = None
        self._interpolate_frames = self.DEFAULT_INTERPOLATE_FRAMES
        self._max_interpolate_distance = self.DEFAULT_MAX_INTERPOLATE_DISTANCE
        self._view_distance=self.DEFAULT_VIEW_DISTANCE
        self._current_direction = self.DEFAULT_DIRECTION
        self.structure.triggers.add_handler('deleted', self._model_deleted_cb)
        if not session_restore:
            get_stepper_mgr(structure.session).register_stepper(self)


    def _block_clipper_spotlights(self):
        from chimerax.clipper import get_all_symmetry_handlers
        sym_handlers = get_all_symmetry_handlers(self.session)
        for sh in sym_handlers:
            sh.triggers._triggers['spotlight moved'].block()

    def _release_clipper_spotlights(self):
        from chimerax.clipper import get_all_symmetry_handlers
        sym_handlers = get_all_symmetry_handlers(self.session)
        for sh in sym_handlers:
            sh.triggers._triggers['spotlight moved'].release()
            if sh.spotlight_mode:
                sh.update()


    def incr_residue(self, direction=None, polymeric_only=True):
        if direction is not None:
            if direction not in ('next', 'prev'):
                from chimerax.core.errors import UserError
                raise UserError('Invalid direction argument! Must be either "next" '
                    'or "prev".')
            self._current_direction = direction
        else:
            direction = self._current_direction
        if direction=='next':
            incr = 1
        else:
            incr = -1
        if polymeric_only:
            residues = self._polymeric_residues()
        else:
            residues = self.structure.residues
        cr = self._current_residue
        if cr is None or cr.deleted:
            next_res = self._first_res(residues, incr)
        else:
            ci = residues.index(cr)
            if (ci == -1 or (ci == len(residues)-1 and incr==1) or
                    (ci==0 and incr==-1)):
                next_res = self._first_res(residues, incr)
            else:
                next_res = residues[ci+incr]
        self._current_residue = next_res
        self._new_camera_position(next_res)
        return next_res

    def _model_deleted_cb(self, *_):
        get_stepper_mgr(self.session).remove_stepper(self)

    def reset_state(self, session):
        self._current_residue=None
        self._interpolate_frames = self.DEFAULT_INTERPOLATE_FRAMES
        self._max_interpolate_distance = self.DEFAULT_MAX_INTERPOLATE_DISTANCE
        self._view_distance=self.DEFAULT_VIEW_DISTANCE
        self._current_direction = self.DEFAULT_DIRECTION

    def incr_chain(self, direction="next"):
        r = self._current_residue
        c = r.chain
        m = self.structure
        i = m.chains.index(c)
        if direction=='next':
            new_i = i+1
            if new_i >= len(m.chains):
                new_i = 0
            new_r = m.chains[new_i].existing_residues[0]
        elif direction=='prev':
            if i == 0:
                new_i = -1
            else:
                new_i = i-1
            new_r = m.chains[new_i].existing_residues[-1]
        else:
            raise TypeError('Direction should be one of "next" or "prev"')
        self.step_to(new_r)

    def next_residue(self, polymeric_only=True):
        return self.incr_residue('next', polymeric_only)

    def previous_residue(self, polymeric_only=True):
        return self.incr_residue('prev', polymeric_only)

    def _go_to_first_residue(self, incr, polymeric_only=True):
        if polymeric_only:
            residues = self._polymeric_residues()
        else:
            residues = self.structure.residues
        r = self._current_residue = self._first_res(residues, incr)
        self._new_camera_position(r)
        return r

    def first_residue(self, polymeric_only=True):
        return self._go_to_first_residue(1, polymeric_only)

    def last_residue(self, polymeric_only=True):
        return self._go_to_first_residue(-1, polymeric_only)

    def step_to(self, residue):
        self._current_residue = residue
        self._new_camera_position(residue)


    def _first_res(self, residues, incr):
        if incr == 1:
            return residues[0]
        return residues[-1]

    def _new_camera_position(self, residue, block_spotlight=True):
        session = self.session
        r = residue
        from chimerax.atomic import Residue, Atoms
        pt = residue.polymer_type
        if pt == Residue.PT_NONE:
            # No preferred orientation
            ref_coords = None
            target_coords = r.atoms.coords
            centroid = target_coords.mean(axis=0)
        elif pt == Residue.PT_AMINO:
            ref_coords = self.peptide_ref_coords
            try:
                target_coords = Atoms([r.find_atom(name) for name in ('N', 'CA', 'C')]).coords
                centroid = target_coords[1]
            except ValueError:
                # Either a key atom is missing, or this is a special residue
                # e.g. NH2
                ref_coords=None
                target_coords = r.atoms.coords
                centroid = target_coords.mean(axis=0)
        elif pt == Residue.PT_NUCLEIC:
            ref_coords = self.nucleic_ref_coords
            try:
                target_coords = Atoms(
                    [r.find_atom(name) for name in ("C2'", "C1'", "O4'")]
                ).coords
                centroid=target_coords[1]
            except ValueError:
                ref_coords=None
                target_coords = r.atoms.coords
                centroid = target_coords.mean(axis=0)

        c = session.main_view.camera
        cp = c.position
        old_cofr = session.main_view.center_of_rotation


        if ref_coords is not None:
            from chimerax.geometry import align_points, Place
            p, rms = align_points(ref_coords, target_coords)
        else:
            from chimerax.geometry import Place
            p = Place(origin=centroid)

        tc = self._camera_ref_pos(self.view_distance)
        np = p*tc
        new_cofr = centroid
        if c.name=='orthographic':
            fw = c.field_width
        else:
            fw = None
        new_fw = self._view_distance*2

        def interpolate_camera(session, f, cp=cp, np=np, oc=old_cofr, nc=new_cofr, fw=fw, nfw=new_fw, vr=self._view_distance, center=np.inverse()*centroid, frames=self._interpolate_frames):
            import numpy
            frac = (f+1)/frames
            v = session.main_view
            c = v.camera
            p = np if f+1==frames else cp.interpolate(np, center, frac=frac)
            cofr = oc+frac*(nc-oc)
            c.position = p
            origin = p.origin()
            vd = c.view_direction()
            dist = numpy.linalg.norm(cofr-origin)
            cp = v.clip_planes
            ncp, fcp = _get_clip_points(session, dist)
            cp.set_clip_position('near', ncp, v)
            cp.set_clip_position('far', fcp, v)
            if c.name=='orthographic':
                c.field_width = fw+frac*(nfw-fw)

        from chimerax.geometry import distance
        if distance(new_cofr, old_cofr) < self._view_distance:
            if block_spotlight:
                self._block_clipper_spotlights()
                from .delayed_reaction import call_after_n_events
                call_after_n_events(self.session.triggers, 'frame drawn', self._interpolate_frames, self._release_clipper_spotlights, [])
            from chimerax.core.commands import motion
            motion.CallForNFrames(interpolate_camera, self._interpolate_frames, session)
        else:
            interpolate_camera(session, 0, frames=1)


    def _polymeric_residues(self, strict=False):
        from chimerax.atomic import Residue
        import numpy
        m = self.structure
        residues = m.residues[m.residues.polymer_types!=Residue.PT_NONE]
        if not len(residues):
            if not strict:
                self.session.logger.warning('You have selected '
                'polymeric_only=True, but this model contains no polymeric '
                'residues. Ignoring the argument. To suppress this warning, use '
                'polymeric_only=False')
                residues = self.structure.residues
            else:
                from chimerax.core.errors import UserError
                raise UserError('You have selected '
                'polymeric_only=True, but this model contains no polymeric '
                'residues. To continue, use polymeric_only=False')
        return residues


    @property
    def field_width(self):
        return self._view_distance

    @property
    def view_distance(self):
        return self._view_distance

    @view_distance.setter
    def view_distance(self, distance):
        self._view_distance = distance

    @property
    def interpolate_frames(self):
        return self._interpolate_frames

    @interpolate_frames.setter
    def interpolate_frames(self, n_frames):
        self._interpolate_frames = n_frames

    @property
    def step_direction(self):
        return self._current_direction

    @step_direction.setter
    def step_direction(self, direction):
        self._current_direction = direction

    @property
    def current_residue(self):
        if self._current_residue is None or self._current_residue.deleted:
            return None
        return self._current_residue

    @property
    def peptide_ref_coords(self):
        '''
        Reference N, CA, C coords for a peptide backbone, used for setting camera
        position and orientation.
        '''
        if not hasattr(self, '_peptide_ref_coords'):
            import numpy
            self._peptide_ref_coords = numpy.array([
                [-0.844, 0.063, -0.438],
                [0.326, 0.668, 0.249],
                [1.680, 0.523, -0.433]
            ])
        return self._peptide_ref_coords

    @property
    def nucleic_ref_coords(self):
        '''
        Reference C2', C1', O4' coords for a nucleic acid residue, used for setting
        camera position and orientation.
        '''
        if not hasattr(self, '_nucleic_ref_coords'):
            import numpy
            self._nucleic_ref_coords = numpy.array([
                [0.822, -0.112, 1.280],
                [0,0,0],
                [0.624, -0.795, -0.985]
            ])
        return self._nucleic_ref_coords

    def _camera_ref_pos(self, distance):
        from chimerax.geometry import Place
        return Place(axes=[[1,0,0], [0,0,1],[0,-1,0]], origin=[0,-distance,0])

    def take_snapshot(self, session, flags):
        print('Taking snapshot of stepper: {}'.format(self.structure.name))
        from . import ISOLDE_STATE_VERSION
        data = {
            'version':  ISOLDE_STATE_VERSION,
            'structure':    self.structure,
            'residue':      self.current_residue,
            'direction':    self._current_direction,
            'interp':       self._interpolate_frames,
            'interp_distance':  self._max_interpolate_distance,
            'view_distance':    self._view_distance,
        }
        return data

    @staticmethod
    def restore_snapshot(session, data):
        stepper = get_stepper(data['structure'], session_restore=True)
        stepper.set_state_from_snapshot(session, data)
        session.triggers.add_handler('end restore session', stepper._end_restore_session_cb)
        return stepper
    
    def _end_restore_session_cb(self, *_):
        '''
        In ISOLDE versions 1.7 and older `ResidueStepper` was a `StateManager` subclass,
        and logic errors were leading to the proliferation of detached `ResidueStepper` 
        instances and eventual tracebacks when sessions were repeatedly saved and restored
        after opening, stepping into and closing models. In 1.8 this was rearranged to have 
        a session-level `ResidueStepperMgr(StateManager)` singleton managing the 
        individual `ResidueStepper(State)` instances. This corrects the original issue, but
        session files restored from 1.7 and earlier come back with spurious entries in 
        `session._state_managers` that must be removed.
        '''
        for tag,inst in self.session.state_managers_by_tag(ResidueStepper).items():
            self.session.remove_state_manager(tag)
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def set_state_from_snapshot(self, session, data):
        self._current_residue = data['residue']
        self._current_direction = data['direction']
        self._interpolate_frames = data['interp']
        self._max_interpolate_distance = data['interp_distance']
        self._view_distance = data['view_distance']

def _get_clip_points(session, dist):
    from chimerax.clipper.mousemodes import ZoomMouseMode
    mm = [b.mode for b in session.ui.mouse_modes.bindings if isinstance(b.mode, ZoomMouseMode)]
    c = session.view.camera
    if len(mm):
        zmm = mm[0]
        o = c.position.origin()
        vd = c.view_direction()
        return (zmm.near_clip_point(o, vd*dist, dist), zmm.far_clip_point(o, vd*dist, dist))
    return (o+vd*dist*0.5, o+vd*dist*1.5)