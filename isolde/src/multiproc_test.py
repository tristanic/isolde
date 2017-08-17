import multiprocessing as mp
from multiprocessing import sharedctypes
from math import ceil

# If using the spawn method, foo needs to be defined somewhere that 
# ChimeraX knows where to find and re-import it. 
def foo(bfactors,start,end,out):
    for i in range(1000):
        c = bfactors.copy()
        c.sort()
    out[start:end] = bfactors
    

def run_multiproc_test(atoms, nproc):
    
    # Have to explicitly re-import foo here
    from chimerax.isolde.multiproc_test import foo
    ctx = mp.get_context('spawn')
    import numpy
    natoms = len(atoms)
    arr = sharedctypes.Array('d', natoms, lock = False)
    ret = numpy.empty(natoms, numpy.float32)
    stride = int(ceil(natoms/nproc))
    proclist = []
    bfactors = atoms.bfactors
    for i in range(nproc):
        start = stride*i
        end = start+stride
        if end > natoms:
            end = natoms
        p = ctx.Process(target=foo, args=(bfactors[start:end],start,end,arr))
        proclist.append(p)
        p.start()
    
    for p in proclist:
        p.join()
    ret[:] = arr
    return ret


import sys
mod = sys.modules['__main__']
mod.__spec__ = None
