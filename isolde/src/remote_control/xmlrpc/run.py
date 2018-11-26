#!/usr/bin/env python

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

# Python 2 test making an XMLRPC request to run a ChimeraX command

if __name__ == '__main__':
    xmlrpc_port = 42184

    from xmlrpc.client import ServerProxy
    s = ServerProxy(uri="http://127.0.0.1:%d/RPC2" % xmlrpc_port)
    print('Available methods: {}'.format(s.system.listMethods()))
    print(s.system.methodHelp('run_command'))
    from sys import argv
    #status = s.run_command(argv[1])
    import os
    arg = argv[1]
    if arg == 'quit':
        status = s.quit()
    else:
        status = s.load_phenix_refine_final_files(os.path.abspath(os.curdir), arg)

    print('Ran "%s", status %s' % (argv[1], status))
