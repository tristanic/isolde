#!/usr/bin/env python

if __name__ == '__main__':
    xmlrpc_port = 42184

    from xmlrpc.client import ServerProxy
    s = ServerProxy(uri="http://127.0.0.1:%d/RPC2" % xmlrpc_port)
    print('Available methods: {}'.format(s.system.listMethods()))
    print(s.system.methodHelp('run_command'))
    from sys import argv
    #status = s.run_command(argv[1])
    import os
    arg = None
    if len(argv) > 1:
        arg = argv[1]
    if arg == 'quit':
        status = s.quit()
    elif arg:
        status = s.load_phenix_refine_final_files(os.path.abspath(os.curdir), arg)

    print('Ran "%s", status %s' % (argv[1], status))
