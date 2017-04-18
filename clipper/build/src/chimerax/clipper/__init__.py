from .main import *
from .crystal import read_mtz
from .clipper_mtz_new import ReflectionDataContainer

def initialize_mouse_modes(session):
    from .mousemodes import ZoomMouseMode, SelectVolumeToContour, ContourSelectedVolume
    z = ZoomMouseMode(session)
    s = SelectVolumeToContour(session)
    v = ContourSelectedVolume(session, s, True)
    session.ui.mouse_modes.bind_mouse_mode('right',[],z)
    session.ui.mouse_modes.bind_mouse_mode('wheel',['control'], s)
    session.ui.mouse_modes.bind_mouse_mode('wheel',[], v) 
