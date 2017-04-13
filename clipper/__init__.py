from .main import *
from .crystal import read_mtz
from .clipper_mtz_new import ReflectionDataContainer

'''
Named triggers to simplify handling of changes on key events (e.g. live
updating of maps/symmetry).
'''
from chimerax.core import triggerset
triggers = triggerset.TriggerSet()
trigger_names = (
    'map box changed',  # Changed shape of box for map viewing
    'map box moved',    # Changed location of box for map viewing
    'atom box changed', # Changed shape of box for showing symmetry atoms
    'atom box moved',   # Changed location of box for showing symmetry atoms
)
for t in trigger_names:
    triggers.add_trigger(t)

def initialize_mouse_modes(session):
    from .mousemodes import ZoomMouseMode, SelectVolumeToContour, ContourSelectedVolume
    z = ZoomMouseMode(session)
    s = SelectVolumeToContour(session)
    v = ContourSelectedVolume(session, s, True)
    session.ui.mouse_modes.bind_mouse_mode('right',[],z)
    session.ui.mouse_modes.bind_mouse_mode('wheel',['control'], s)
    session.ui.mouse_modes.bind_mouse_mode('wheel',[], v) 
