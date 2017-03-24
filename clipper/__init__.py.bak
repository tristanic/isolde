from chimerax.core.toolshed import BundleAPI


class _MyAPI(BundleAPI):

    @staticmethod
    def get_class(class_name):
        # 'get_class' is called by session code to get class from bundle that
        # was saved in a session
        if class_name == 'Clipper_ToolInstance':
            from . import tool
            return tool.Clipper_ToolInstance
        return None

    @staticmethod
    def start_tool(session, tool_name):
        # 'start_tool' is called to start an instance of the tool
        from .tool import Clipper_ToolInstance
        from chimerax.core import tools
        return tools.get_singleton(session, Clipper_ToolInstance, 'Clipper', create=True)
        

    @staticmethod
    def register_command(command_name):
        # 'register_command' is lazily called when the command is referenced
        #from chimerax.core.commands import register
        pass

bundle_api = _MyAPI()
