import multiprocessing as mp
from multiprocessing import sharedctypes, Pool, spawn
import ctypes
from math import ceil
import numpy
from numpy import ndarray
from .threading import SharedNumpyArray

try:
    from chimerax import app_bin_dir
    import os
    spawn.set_executable(os.path.join(app_bin_dir, 'python3.6'))
except:
    # We're not in ChimeraX any more, Toto!
    pass

def error_callback(e):
    print(e)

 
#~ class SharedNumpyArray(ndarray):
    #~ '''
    #~ multiprocessing.Array types are thread-safe by default, but are 
    #~ horribly inefficient in getting/setting data. If you want speed you
    #~ need to create a Numpy array pointing to the same shared memory,
    #~ but this circumvents the automatic acquire/release behaviour. To 
    #~ provide thread-safe behaviour would therefore require carrying through
    #~ both the original Array (for the lock) and the derived Numpy array.
    #~ This class is an attempt to get the best of both worlds: it behaves
    #~ just like a normal Numpy array, but carries through the Array lock
    #~ object and its methods. To use it:
    
    #~ (In master process):
    #~ import multiprocessing as mp
    #~ mp_array = mp.Array(type, data_or_init)
    #~ shared_numpy = SharedNumpyArray(mp_array)
    
      #~ Pass shared_numpy to the thread Pool init function or to the thread
      #~ itself if creating threads on the fly.
    
    #~ (In each thread):
      #~ If thread safety is not required (that is, different threads don't
      #~ attempt to read and write to the same index), then just use it like
      #~ any other array. If thread safety *is* required:
    
    #~ with shared_numpy.get_lock():
        #~ do_something(shared_numpy)
    #~ '''
    #~ def __new__(cls, mp_array):
        #~ if mp_array is None:
            #~ raise TypeError('Please provide a multiprocessing.Array object\
                             #~ with a thread lock!')
        #~ obj = numpy.frombuffer(mp_array.get_obj(), type(mp_array[0])).view(cls)
        #~ obj._mparray = mp_array
        #~ obj.get_lock = mp_array.get_lock
        #~ obj.acquire = mp_array.acquire
        #~ obj.release = mp_array.release
        #~ return obj

    #~ def __array_finalize__(self, obj):
        #~ if obj is None: 
            #~ return
        
        #~ self._mparray = getattr(obj, '_mparray', None)
        #~ self.get_lock = getattr(obj, 'get_lock', None)
        #~ self.acquire = getattr(obj, 'acquire', None)
        #~ self.release = getattr(obj, 'release', None)
        
        
        

def _pool_init(ret, sca, eca, n_points, n_coords):
    '''
    Creates and sets up the shared variables needed by the threads. Only
    runs within the threads themselves, so the global variables aren't 
    too evil. Not sure if there's a way around using them.
    '''
    global shared_arr
    global start_c
    global end_c
    shared_arr = ret
    start_c = sca
    end_c = eca
    
def _interpolate_worker(interp_fracs,frames, proc_id):
    global shared_arr
    from chimerax.core.geometry import interpolate_points
    import numpy
    global start_c
    global end_c
    target = numpy.empty((len(frames), len(start_c),3), numpy.double)
    
    for (frac, frame) in zip(interp_fracs, frames):
        shared_arr[frame] = interpolate_points(start_c, end_c, frac)
    return True
    

def run_multiproc_test(start_coords, end_coords, num_interpolation_points, nproc, ncycles=2, context = 'spawn'):
    from time import time, sleep
    start_time = time()    
    ctx = mp.get_context(context)
    import numpy
    n_coords = len(start_coords)
    assert n_coords == len(end_coords)
    n_points = num_interpolation_points
    c_arr = mp.Array(ctypes.c_double, n_points*n_coords*3)
    ret = SharedNumpyArray(c_arr).reshape((n_points, n_coords, 3))
    start_coord_array = mp.Array(ctypes.c_double, n_coords*3)
    end_coord_array = mp.Array(ctypes.c_double, n_coords*3)
    sca = SharedNumpyArray(start_coord_array).reshape((n_coords,3))
    sca[:] = start_coords
    eca = SharedNumpyArray(end_coord_array).reshape((n_coords,3))
    eca[:] = end_coords
    frames = numpy.array(range(num_interpolation_points),dtype='int')
    fracs = frames / n_points
    stride = int(ceil(n_points/nproc))
    with Pool(processes = nproc, initializer = _pool_init, 
              initargs = (ret, sca, eca, n_points, n_coords)) as p:
        for cycle in range(ncycles):
            results = []
            for i in range(nproc):
                start = stride*i
                end = start+stride
                if end > n_points:
                    end = n_points
                results.append(p.apply_async(_interpolate_worker,
                    args=(fracs[start:end],frames[start:end], i),
                    error_callback=error_callback))
            count = 0
            while True:
                done = True
                for result in results:
                    if not result.ready():
                        done = False
                if done:
                    break
                count += 1
                if count > 50000:
                    print('Timeout!')
                    break
                sleep(1e-4)
            print('Cycle {} took {} seconds'.format(cycle, time() - start_time))
            start_time = time()
            
    
    
    print('Finishing took {} seconds'.format(time()-start_time))
    return ret

