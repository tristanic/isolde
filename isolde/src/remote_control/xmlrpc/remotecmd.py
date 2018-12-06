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

def register_remote_control_command(command_name, logger):

    from chimerax.core.commands import CmdDesc, register, BoolArg, StringArg, IntArg, FloatArg
    desc = CmdDesc(
        required = [('enable', BoolArg)],
        keyword = [('address', StringArg),
                   ('port', IntArg),
                   ('timeout', FloatArg)],
        synopsis = 'Allow other processes to send XMLRPC commands to ChimeraX')
    register(command_name, desc, remote_control, logger=logger)
    register('remotecontrol_xmlrpc_on', CmdDesc(synopsis='start xmlrpc server'), remote_control_xmlrpc_on, logger=logger)

def remote_control(session, enable, address = '127.0.0.1', port = 42184, timeout = 0.01):
    '''
    Start XMLRPC server so Phenix can remotely control ISOLDE.
    '''
    isolde = getattr(session, 'isolde', None)
    if isolde is None:
        from chimerax.isolde import bundle_api
        bundle_api.start_tool(session, 'ISOLDE')
        isolde = session.isolde
    s = getattr(isolde, '_xmlrpc_server', None)
    if enable:
        if s is None:
            s = Isolde_XMLRPC_Server(isolde, address, port, timeout)
            isolde._xmlrpc_server = s
        return s
    else:
        if s:
            del isolde._xmlrpc_server
            s.close()

def remote_control_xmlrpc_on(session):
    remote_control(session, True)
#
# XML remote procedure call server run by ChimeraX to allow other apps such as
# the Phenix x-ray model building program to use ChimeraX for visualizing results.
#


class Isolde_XMLRPC_Server:
    REGISTERED_METHODS=(
        'run_command',
        'spotlight_mode',
        'is_alive',
        'quit',
        'close_all',
        'close_maps',
        'show_start_model',
        'close_tmp_model',
        'clear_refinement',
        'load_phenix_refine_temp_files',
        'load_phenix_refine_final_files'
    )
    def __init__(self, isolde, address='127.0.0.1', port=42184, timeout=0.01):
        self.isolde = isolde
        self.session = isolde.session
        from queue import Queue
        cq = self._cmd_queue = Queue()
        rq = self._result_queue = Queue()
        st = self._server_thread = _Isolde_XMLRPC_Server_Thread(
            cq, rq, address, port, timeout
        )

        self._current_model = None
        self._start_model = None
        self.settings = { "map_color" : [0.0, 0.5, 1.0, 1.0],
                          "diff_map_colors" : ([0.0, 1.0, 0.0, 1.0],
                                                [1.0, 0.0, 0.0, 1.0]),
                          "n_map_color" : [0.25, 0.0, 1.0, 1.0],
                          "n_diff_map_colors" : ([0.5, 1.0, 0.0, 1.0],
                                                [1.0, 0.5, 0.0, 1.0]),
                          "anom_map_color" : [1.0, 1.0, 0.0, 1.0],
                          "iso_diff_map_colors" : ([0.5, 0.0, 1.0, 1.0],
                                                  [1.0, 0,0, 0.5, 1.0]), }

        self._cmd_check_handler = self.session.triggers.add_handler(
            'new frame', self._thread_cmd_cb
        )
        self._app_quit_handler = self.session.triggers.add_handler(
            'app quit', self.shutdown
        )
        for f_name in self.REGISTERED_METHODS:
            f = getattr(self, f_name)
            st.register_function(f)
        # st.register_function(self.run_command)
        st.start()

    def run_command(self, command):
        '''
        Run an arbitrary ChimeraX command (that is, any command that would be
        valid if typed on the ChimeraX command line). The argument should be
        a single string.
        '''
        from chimerax.core.commands import run
        try:
            run(self.session, command)
        except Exceptions as e:
            self.session.info(str(e))
            return False
        return True

    def is_alive(self):
        return True

    def spotlight_mode(self):
        '''
        Force a return to "spotlight" viewing mode.
        '''
        from chimerax.clipper.symmetry import get_symmetry_handler
        sh = get_symmetry_handler(self._current_model)
        sh.spotlight_mode=True

    def set_default_atom_coloring(self):
        self.run_command('color bychain')
        self.run_command('color byhetero')

    def quit(self):
        '''
        Shut down the ChimeraX session. No arguments.
        '''
        # Give the result to the thread immediately to avoid timeouts
        self._result_queue.put(True)
        self.run_command('exit')
        return True

    def close_all(self):
        '''
        Close all currently open models and maps. No arguments.
        '''
        self.run_command('close')
        self._current_model = None
        self._start_model = None
        return True

    def clear_refinement(self):
        '''
        Close the current refinement result and its maps. No arguments.
        '''
        if self._current_model is not None and not self._current_model.deleted:
            self.session.models.close([self._current_model.parent])
            self._current_model = None
        return True

    def close_tmp_model(self):
        '''
        Synonym for clear_refinement(). No arguments.
        '''
        return self.clear_refinement()

    def show_start_model(self, file_name):
        '''
        Load the starting model from file. Takes the filename (including full
        path) as an argument.
        '''
        if self._start_model is not None and not self._start_model.deleted:
            self.session.models.close([self._start_model])
        from chimerax.core.commands import open as cxopen
        m = self._start_model = cxopen.open(self.session, file_name)[0]
        m.bonds.radii=0.05
        from chimerax.clipper.symmetry import get_symmetry_handler
        get_symmetry_handler(m)
        self.set_default_atom_coloring()
        return True

    def _load_model_and_mtz(self, pdb_file, mtz_file):
        import os
        if (not os.path.isfile(pdb_file)) or (not os.path.isfile(mtz_file)):
            raise RuntimeError('One or more output files not found: \n'
                               '{}\n{}'.format(pdb_file, mtz_file))
        self.clear_refinement()
        from chimerax.core.commands import open as cxopen
        m = self._current_model = cxopen.open(self.session, pdb_file)[0]
        self.isolde.add_xtal_data(mtz_file, model=m)
        self.set_default_atom_coloring()
        self.isolde.change_selected_model(m)


    def load_phenix_refine_temp_files(self, tmp_dir, run_name):
        '''
        Load the most recent interim output from phenix.refine. Will
        automatically close the last output. Takes two arguments: the path to
        the PHENIX temp directory, and the base name for this refinement run.
        '''
        import os
        if None in [tmp_dir, run_name] or not os.path.isdir(tmp_dir):
            raise RuntimeError('Invalid arguments or nonexistent directory!')
        pdb_tmp = os.path.join(tmp_dir, 'tmp.refine.pdb')
        map_tmp = os.path.join(tmp_dir, 'tmp.refine_maps.mtz')
        self._load_model_and_mtz(pdb_tmp, map_tmp)
        return True

    def _load_phenix_refine_final_files(self, output_dir, file_base,
            joint_xray_neutron=False):
        import os
        if None in [output_dir, file_base] or not os.path.isdir(output_dir):
            raise RuntimeError('Invalid arguments or nonexistent directory!')
        output_base = os.path.join(output_dir, file_base)

        pdb_file = os.path.join(output_base+'.pdb')
        map_file_old = output_base+'_map_coeffs.mtz'
        map_file_new = output_base+'.mtz'
        map_file = None
        if os.path.isfile(map_file_new):
            map_file = map_file_new
        elif os.path.isfile(map_file_old):
            map_file=map_file_old
        if map_file is not None:
            self._load_model_and_mtz(pdb_file, map_file)
        else:
            raise RuntimeError('Could not find MTZ file for run name {}'.format(
                file_base
            ))

    def load_phenix_refine_final_files(self, output_dir, file_base):
        '''
        Load the final phenix.refine results for a standard x-ray refinement.
        Will automatically close the temporary files.
        '''
        self._load_phenix_refine_final_files(output_dir, file_base,
            joint_xray_neutron=False)
        return True


    def load_phenix_refine_final_xn_files(self, output_dir, file_base):
        '''
        Load the final phenix.refine results for a joint x-ray/neutron
        refinement. Will automatically close the temporary files.
        '''
        self._load_phenix_refine_final_files(output_dir, file_base,
            joint_xray_neutron=True)
        m = self._current_model
        from chimerax.clipper.symmetry import get_symmetry_handler
        sh = get_symmetry_handler(m)
        xmapset = sh.map_mgr.xmapsets[0]
        return True

    def recenter_and_zoom(self, x, y, z):
        import numpy
        self.spotlight_mode()
        from chimerax.isolde.view import focus_on_coord
        focus_on_coord(self.session, numpy.array([x,y,z]), radius=5.0)
        return True


    def close_maps(self):
        if self._current_model is None or self._current_model.deleted:
            raise RuntimeError('No PHENIX model/maps currently open!')
        from chimerax.clipper.symmetry import get_symmetry_handler
        sh = get_symmetry_handler(self._current_model)
        if len(sh.map_mgr.xmapsets):
            self.session.models.close(sh.map_mgr.xmapsets)
        return True

    # NON-EXPOSED METHODS
    def shutdown(self, *_):
        self._server_thread.server.shutdown()
        self._server_thread.join()

    # CALLBACKS

    def _thread_cmd_cb(self, *_):
        cq = self._cmd_queue
        if not cq.empty():
            func_name, args = cq.get()
            f = getattr(self, func_name)
            try:
                result = f(*args)
            except Exception as err:
                self.session.logger.warning('ISOLDE XMLRPC error in {}: \n'.format(f.__name__)
                    + str(err))
                result = False
            self._result_queue.put(result)


def _rpc_method_factory(cmd_queue, result_queue, func):
    import functools
    @functools.wraps(func)
    def f(*args, **kwargs):
        if not result_queue.empty():
            result_queue.clear()
        cmd_queue.put((func.__name__, args))
        try:
            result = result_queue.get(timeout=10)
        except:
            result = False
        if result is None:
            result = False
        return result
    return f


from threading import Thread
from xmlrpc.server import SimpleXMLRPCServer
class _Isolde_XMLRPC_Server_Thread(Thread):

    def __init__ (self, cmd_queue, result_queue, address, port, timeout):
        super().__init__()

        self.cmd_queue = cmd_queue
        self.result_queue = result_queue
        # start XML-RPC server
        server = self._server = SimpleXMLRPCServer((address, port), logRequests=0)
        server.socket.settimeout(timeout)
        server.register_introspection_functions()

    def run(self):
        self._server.serve_forever()


    def register_function(self, function):
        f = _rpc_method_factory(self.cmd_queue, self.result_queue, function)
        self._server.register_function(f, name=function.__name__)

    @property
    def server(self):
        return self._server
