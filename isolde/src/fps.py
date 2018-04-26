# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



def fps_start(session):
    singleton(session).start_showing()
from chimerax.core.commands import CmdDesc
fps_start_desc = CmdDesc()

def fps_stop(session):
    singleton(session).stop_showing()
from chimerax.core.commands import CmdDesc
fps_stop_desc = CmdDesc()

_fps_tracker = None
def singleton(session):
    global _fps_tracker
    if _fps_tracker is None:
        _fps_tracker = Track_FPS(session)
    return _fps_tracker


class Track_FPS():

    def __init__(self, session):
        self.sess_triggers = session.triggers
        self.log = session.logger.status
        self.showing = False
        self.last_time = None
        self.text = None
        self.counter = 0
        self.num_frames_to_average = 4

    def start_showing(self):
        if self.showing:
            return True
        self.handler = self.sess_triggers.add_handler('new frame',
            self.update_fps)
        from time import time
        self.last_time = time()
        self.showing = True

    def stop_showing(self):
        self.sess_triggers.remove_handler(self.handler)
        self.log('')
        self.showing = False

    def update_fps(self, *_):
        self.counter = (self.counter + 1) % self.num_frames_to_average
        if self.counter != 0:
            return
        from time import time
        current_time = time()
        time_between_frames = current_time - self.last_time
        if time_between_frames != 0:
            fps = 1 / time_between_frames * self.num_frames_to_average
        else:
            fps = 0
        fps_str = '{:0.1f}'.format(fps)
        self.last_time = current_time
        display_text = 'Framerate ' + fps_str + ' frames/sec'
        self.log(display_text)
