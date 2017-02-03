# vim: set expandtab shiftwidth=4 softtabstop=4:

# ToolUI should inherit from ToolInstance if they will be
# registered with the tool state manager.
#
# ToolUI classes may also override
#   "delete" - called to clean up before instance is deleted
#
from chimerax.core.tools import ToolInstance


class Clipper_ToolInstance(ToolInstance):

    SESSION_ENDURING = False
    # if SESSION_ENDURING is True, tool instance not deleted at session closure

    def __init__(self, session, tool_name):
        ToolInstance.__init__(self, session, tool_name)
        self.display_name = "Clipper"
        
        from . import main
        self.main = main
        session.clipper = self.main
        session.clipper.session = session
        
    
