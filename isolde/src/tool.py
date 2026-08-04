# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 05-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



# ToolUI should inherit from ToolInstance if they will be
# registered with the tool state manager.
#
# ToolUI classes may also override
#   "delete" - called to clean up before instance is deleted
#
from chimerax.core.tools import ToolInstance
from Qt.QtCore import Qt 

class ISOLDE_ToolUI(ToolInstance):
    SESSION_ENDURING = True

    # Target of the "Help" entry in the tool window's right-click context menu
    # (ToolInstance.help defaults to None, which makes that entry read "No help
    # available").  This points at the GUI's own doc page, which is narrower than
    # the docs root the GUI's Help button opens via
    # Isolde.show_master_help_in_browser.  Built by make_docs.bat into
    # src/docs/user/, which ChimeraX serves as help:user/.
    help = 'help:user/tools/ISOLDE.html'

    def __init__(self, session, tool_name, show_splash=True):
        super().__init__(session, tool_name)
        from .isolde import Isolde
        isolde = Isolde(session)

        from chimerax.toolbar.tool import get_toolbar_singleton
        tb = get_toolbar_singleton(session, create=False)
        if tb is not None:
            tb.ttb.show_tab('ISOLDE')

        self.display_name='ISOLDE'
        if show_splash:
            self._show_splash()
        self.session.triggers.add_handler('new frame', self._launch_main_gui)
    
    def _launch_main_gui(self, *_):
        from .ui.main_win import IsoldeMainWin
        tw = self.tool_window = IsoldeMainWin(self)
        from chimerax.log.tool import Log
        # log_instances = self.session.tools.find_by_class(Log)
        # if len(log_instances):
        #     placement = log_instances[0].tool_window
        #     # Make sure it appears as the front-most tab
        #     self.session.triggers.add_handler('new frame', self._raise_tab_to_front)
        # else:
        #     placement='side'
        tw.manage(placement=None, allowed_areas=Qt.LeftDockWidgetArea|Qt.RightDockWidgetArea)
        tw.shrink_to_fit()
        isolde = self.session.isolde
        isolde.triggers.activate_trigger(isolde.GUI_STARTED, None)
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER
    
    def _raise_tab_to_front(self, *_):
        tw = self.tool_window
        tw.ui_area.parent().parent().raise_()
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER



    def _show_splash(self):
        from Qt.QtGui import QPixmap
        from Qt.QtWidgets import QSplashScreen
        from Qt.QtCore import Qt
        import os
        session = self.session
        root_dir = os.path.dirname(os.path.abspath(__file__))
        splash_pix = QPixmap(os.path.join(
            root_dir,'resources/isolde_splash_screen.jpg'))
        splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
        splash.setMask(splash_pix.mask())
        splash.show()
        # Register the splash so any blocking dialog raised during startup
        # (e.g. the altloc or unbonded-disulfide warnings) can dismiss it first
        # rather than be hidden behind this always-on-top window.
        from .dialog import register_splash
        register_splash(splash)
        from time import time
        start_time = [time()]
        def _splash_remove_cb(trigger_name, data, start_time=start_time, min_time=2):
            from time import time
            elapsed_time = time()-start_time[0]
            if elapsed_time > min_time:
                start_time[0] = time()
                session.triggers.add_handler('new frame', _splash_fade_cb)
                from chimerax.core.triggerset import DEREGISTER
                return DEREGISTER
        def _splash_fade_cb(trigger_name, data, splash=splash, start_time=start_time, fade_time=0.25):
            from time import time
            et = time()-start_time[0]
            opacity = 1-et/fade_time
            if opacity <= 0:
                splash.close()
                from .dialog import register_splash
                register_splash(None)
                from chimerax.core.triggerset import DEREGISTER
                return DEREGISTER
            splash.setWindowOpacity(opacity)
        session.triggers.add_handler('new frame', _splash_remove_cb)



    def delete(self):
        self.session.isolde._on_close()
        self.tool_window.cleanup()
        super().delete()
