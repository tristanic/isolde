import os, sys, glob
import numpy
import ctypes
from chimerax.core.state import State
from chimerax.core.atomic import molc
from chimerax.core.atomic.molc import CFunctions, string, cptr, pyobject, \
    set_c_pointer, pointer, size_t
# object map lookups
from chimerax.core.atomic.molobject import _atoms, \
                _atom_pair, _atom_or_none, _bonds, _chain, _element, \
                _pseudobonds, _residue, _residues, _rings, _non_null_residues, \
                _residue_or_none, _residues_or_nones, _residues_or_nones, \
                _chains, _atomic_structure, _pseudobond_group, \
                _pseudobond_group_map

from numpy import int8, uint8, int32, uint32, float64, float32, byte, bool as npy_bool

libdir = os.path.dirname(os.path.abspath(__file__))
libfile = glob.glob(os.path.join(libdir, '..', 'openmm.cpython*'))[0]

_c_functions = CFunctions(os.path.splitext(libfile)[0])
c_property = _c_functions.c_property
cvec_property = _c_functions.cvec_property
c_function = _c_functions.c_function
c_array_function = _c_functions.c_array_function

class OpenMM_Thread_Handler:
    '''
    A lightweight wrapper class for an OpenMM simulation `Context`, which pushes
    the task of stepping the simulation forward off to a C++ thread so that
    Python performance is not interrupted. Each call to `step()` or `minimize()`
    will create a short-lived thread to run the desired number of steps or
    a round of minimization, respectively. The status of the thread can be
    checked with thread_finished, while the initial and final coordinates can
    be retrieved with `last_coords` and `coords` respectively. If instability
    (overly fast-moving atoms) is detected, the thread will terminate early and
    the `unstable` property will set to True. In such cases it is advisable to
    run one or more minimization rounds. When minimization converges to with
    tolerance, `unstable` will be reset to False.

    Use with care! While nothing prevents you from using the standard single-
    thread OpenMM API alongside this one, it is up to you to ensure that
    no threads are running before making any calls that affect the simulation.
    Additionally, the `OpenMM_Thread_Handler` object *must* be destroyed before
    the Context it is attached to.
    '''
    def __init__(self, context, c_pointer=None):
        cname = type(self).__name__.lower()
        if c_pointer is None:
            new_func = cname + '_new'
            int_ptr = int(context.this)
            f = c_function(new_func, args=(ctypes.c_void_p,), ret=ctypes.c_void_p)
            c_pointer = f(int_ptr)
        set_c_pointer(self, c_pointer)
        f = c_function('set_'+cname+'_py_instance', args=(ctypes.c_void_p, ctypes.py_object))
        f(self._c_pointer, self)
        self.context = context

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def delete(self):
        c_function('openmm_thread_handler_delete', args=(ctypes.c_void_p,))(self._c_pointer)

    def step(self, steps):
        f = c_function('openmm_thread_handler_step',
            args=(ctypes.c_void_p, ctypes.c_size_t))
        f(self._c_pointer, steps)

    def minimize(self):
        '''
        Run an energy minimization on the coordinates. If the minimisation
        converges to within tolerance, unstable will be set to False.
        Don't forget to run reinitialize_velocities() before continuing
        equilibration!
        '''
        f = c_function('openmm_thread_handler_minimize',
            args = (ctypes.c_void_p,))
        f(self._c_pointer)

    def reinitialize_velocities(self):
        if not self.thread_finished:
            raise RuntimeError('Wait for the thread to finish before reinitialising velocities!')
        c = self.context
        i = c.getIntegrator()
        c.setVelocitiesToTemperature(i.getTemperature())

    @property
    def natoms(self):
        f = c_function('openmm_thread_handler_num_atoms',
            args=(ctypes.c_void_p,),
            ret = ctypes.c_size_t)
        return f(self._c_pointer)

    def thread_finished(self):
        f = c_function('openmm_thread_handler_thread_finished',
            args=(ctypes.c_void_p,),
            ret=npy_bool)
        return f(self._c_pointer)

    def finalize_thread(self):
        f = c_function('openmm_thread_handler_finalize_thread',
            args=(ctypes.c_void_p,))
        f(self._c_pointer)

    def unstable(self):
        f = c_function('openmm_thread_handler_unstable',
            args=(ctypes.c_void_p,),
            ret = ctypes.c_bool)
        return f(self._c_pointer)

    @property
    def unstable_atoms(self):
        f = c_function('openmm_thread_handler_unstable_atoms',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p))
        n = self.natoms
        ret = numpy.zeros(n, numpy.bool)
        f(self._c_pointer, n, pointer(ret))
        return ret

    @property
    def last_coords(self):
        f = c_function('openmm_thread_handler_last_coords',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p))
        n = self.natoms
        coords = numpy.empty((n,3), float64)
        f(self._c_pointer, n, pointer(coords))
        return coords

    @property
    def coords(self):
        f = c_function('openmm_thread_handler_current_coords',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p))
        n = self.natoms
        coords = numpy.empty((n,3), float64)
        f(self._c_pointer, n, pointer(coords))
        return coords

    @coords.setter
    def coords(self, coords):
        f = c_function('set_openmm_thread_handler_current_coords',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p))
        n = self.natoms
        f(self._c_pointer, n, pointer(coords))
        self.reinitialize_velocities()

    def _get_min_thread_period(self):
        '''Throttle the simulation to a minimum time period per loop (in ms)'''
        f = c_function('openmm_thread_handler_min_thread_period',
            args=(ctypes.c_void_p,),
            ret=ctypes.c_double)
        return f(self._c_pointer)

    def _set_min_thread_period(self, period):
        f = c_function('set_openmm_thread_handler_min_thread_period',
            args=(ctypes.c_void_p, ctypes.c_double))
        f(self._c_pointer, period)

    min_thread_period = property(_get_min_thread_period, _set_min_thread_period)
