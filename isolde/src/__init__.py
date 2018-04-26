# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



from chimerax.core.toolshed import BundleAPI
from . import geometry


class _MyAPI(BundleAPI):

    @staticmethod
    def get_class(class_name):
        # 'get_class' is called by session code to get class from bundle that
        # was saved in a session
        if class_name == 'ISOLDE_ToolUI':
            from . import tool
            return tool.ISOLDE_ToolUI
        return None

    @staticmethod
    def start_tool(session, tool_name):
        # 'start_tool' is called to start an instance of the tool
        from .tool import ISOLDE_ToolUI
        from chimerax.core import tools
        return tools.get_singleton(session, ISOLDE_ToolUI, 'ISOLDE', create=True)


    @staticmethod
    def register_command(command_name, logger):
        # 'register_command' is lazily called when the command is referenced
        from . import fps
        from chimerax.core.commands import register
        if command_name == 'fps':
            register(command_name + " start", fps.fps_start_desc, fps.fps_start, logger)
            register(command_name + " stop", fps.fps_stop_desc, fps.fps_stop, logger)
        elif command_name == 'isolde':
            from . import cmd
            cmd.register_isolde()

bundle_api = _MyAPI()
