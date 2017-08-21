import multiprocessing as mp
from multiprocessing import sharedctypes, Pool, spawn
import ctypes
from math import ceil

try:
    from chimerax import app_bin_dir
    import os
    spawn.set_executable(os.path.join(app_bin_dir, 'python3.6'))
except:
    # We're not in ChimeraX any more, Toto!
    pass

def error_callback(e):
    print(e)

def _pool_init(c_arr, start_coord_array, end_coord_array, n_points, n_coords):
    '''
    Creates and sets up the shared variables needed by the threads. Only
    runs within the threads themselves, so the global variables aren't 
    too evil. Not sure if there's a way around using them.
    '''
    import numpy
    global shared_arr
    global start_c
    global end_c
    shared_arr = numpy.frombuffer(c_arr.get_obj()).reshape((n_points, n_coords, 3))
    start_c = numpy.frombuffer(start_coord_array.get_obj()).reshape((n_coords, 3))
    end_c = numpy.frombuffer(end_coord_array.get_obj()).reshape((n_coords, 3))
    
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
    start_coord_array = mp.Array(ctypes.c_double, n_coords*3)
    end_coord_array = mp.Array(ctypes.c_double, n_coords*3)
    sca = numpy.frombuffer(start_coord_array.get_obj()).reshape((n_coords,3))
    sca[:] = start_coords
    eca = numpy.frombuffer(end_coord_array.get_obj()).reshape((n_coords,3))
    eca[:] = end_coords
    frames = numpy.array(range(num_interpolation_points),dtype='int')
    fracs = frames / n_points
    stride = int(ceil(n_points/nproc))
    with Pool(processes = nproc, initializer = _pool_init, 
              initargs = (c_arr,start_coord_array,end_coord_array, n_points, n_coords)) as p:
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
            
    
    
    ret = numpy.frombuffer(c_arr.get_obj()).reshape((n_points, n_coords, 3))
    print('Finishing took {} seconds'.format(time()-start_time))
    return ret

