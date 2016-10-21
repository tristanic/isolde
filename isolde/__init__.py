# vim: set expandtab ts=4 sw=4:

# === UCSF ChimeraX Copyright ===
# Copyright 2016 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  For details see:
# http://www.rbvi.ucsf.edu/chimerax/docs/licensing.html
# This notice must be embedded in or attached to all copies,
# including partial copies, of the software or any revisions
# or derivations thereof.
# === UCSF ChimeraX Copyright ===

from chimerax.core import commands

commands.camera.camera(session,type='ortho')

from fps import Track_FPS
session.fpstracker = Track_FPS(session)
session.fpstracker.start_showing()



from isolde import isolde
session.isolde = isolde.Isolde(session)
session.ui.main_window.add_custom_menu_entry('Plugins','ISOLDE',session.isolde.start_gui)
#session.isolde.start_gui()
#session.isolde.mainwin.show()


