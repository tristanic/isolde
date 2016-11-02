# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.toolshed import BundleAPI


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
        return ISOLDE_ToolUI(session, tool_name)

    @staticmethod
    def register_command(command_name):
        # 'register_command' is lazily called when the command is referenced
        from . import cmd
        from chimerax.core.commands import register
        register(command_name + " start", cmd.fps_start_desc, cmd.fps_start)
        register(command_name + " stop", cmd.fps_stop_desc, cmd.fps_stop)

bundle_api = _MyAPI()
