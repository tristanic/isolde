# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 11-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def register_remote_commands(logger):
    from .xmlrpc.remotecmd import register_isolde_xmlrpc_server
    register_isolde_xmlrpc_server(logger)
    from.rest_server.cmd import register_isolde_rest_server
    register_isolde_rest_server(logger)
