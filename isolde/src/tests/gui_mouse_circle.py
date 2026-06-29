# @Author: Tristan Croll
# @Date:   19-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 19-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
GUI-ONLY demo: move the mouse cursor in a circle for ~5 seconds.

This is NOT a pytest test -- it needs a Qt event loop, so it cannot run under
``--nogui``. Run it inside a normal GUI ChimeraX session:

    runscript /path/to/isolde/src/tests/gui_mouse_circle.py

It just drives the physical cursor (QCursor.setPos) around a circle centred on
the ChimeraX main window. No model, simulation, or mouse mode involved.

NOTE: written but NOT executed by the author (this environment is headless), so
it has only been syntax-checked.
'''

import math

DURATION_S = 5.0
RADIUS_FRAC = 0.15       # circle radius as a fraction of the smaller window dim
TICK_MS = 16             # ~60 cursor updates per second
PERIOD_S = 1.25          # seconds per revolution -> ~4 revolutions in 5 s


class MouseCircleDriver:
    def __init__(self, session):
        self.session = session
        self.log = session.logger

        from Qt.QtGui import QGuiApplication
        # Centre the circle on the ChimeraX main window if we can, else on the
        # primary screen.
        cx = cy = radius = None
        try:
            mw = session.ui.main_window
            from Qt.QtCore import QPoint
            w, h = mw.width(), mw.height()
            centre = mw.mapToGlobal(QPoint(w // 2, h // 2))
            cx, cy = centre.x(), centre.y()
            radius = RADIUS_FRAC * min(w, h)
        except Exception:
            geo = QGuiApplication.primaryScreen().geometry()
            cx, cy = geo.center().x(), geo.center().y()
            radius = RADIUS_FRAC * min(geo.width(), geo.height())

        self.cx, self.cy, self.radius = cx, cy, radius
        self.elapsed = 0.0
        self.steps = 0
        self.log.info('Mouse-circle: orbiting cursor for {:.0f} s '
                      '(centre {},{}; radius {:.0f} px).'.format(
                          DURATION_S, int(cx), int(cy), radius))
        self._start_timer()

    def _start_timer(self):
        from Qt.QtCore import QTimer
        # Parent to the main window so Qt keeps the timer alive after the
        # runscript scope exits (otherwise it is garbage-collected before it
        # ever fires).
        parent = getattr(self.session.ui, 'main_window', None)
        self._timer = QTimer(parent)
        self._timer.timeout.connect(self._tick)
        self._timer.start(TICK_MS)

    def _tick(self):
        self.elapsed += TICK_MS / 1000.0
        if self.elapsed >= DURATION_S:
            self._timer.stop()
            self.log.info('Mouse-circle done: {} cursor moves over {:.1f} s.'.format(
                self.steps, self.elapsed))
            # Release the reference we stashed on the session.
            if getattr(self.session, '_mouse_circle_driver', None) is self:
                self.session._mouse_circle_driver = None
            return
        from Qt.QtCore import QPoint
        from Qt.QtGui import QCursor
        theta = 2 * math.pi * (self.elapsed / PERIOD_S)
        x = self.cx + self.radius * math.cos(theta)
        y = self.cy + self.radius * math.sin(theta)
        QCursor.setPos(QPoint(int(x), int(y)))
        self.steps += 1


def run_mouse_circle(session):
    if not session.ui.is_gui:
        session.logger.warning(
            'gui_mouse_circle requires GUI mode; it cannot run under --nogui.')
        return None
    # Stash on the session so the driver (and its timer) survive the runscript
    # scope; without a persistent reference it is GC'd before the timer fires.
    driver = MouseCircleDriver(session)
    session._mouse_circle_driver = driver
    return driver


# When launched via "runscript", ChimeraX injects `session` into globals.
if 'session' in dir():
    run_mouse_circle(session)        # noqa: F821
