from numpy import ndarray, frombuffer
import multiprocessing as mp
from multiprocessing import sharedctypes
import ctypes


def TypedMPArray(typecode_or_type, size_or_initializer, *, lock=True, ctx=None):
    obj = mp.Array(typecode_or_type, size_or_initializer)
    
    obj.dtype = mp.sharedctypes.typecode_to_type.get(typecode_or_type, typecode_or_type)
    
    return obj
    
    


class SharedNumpyArray(ndarray):
    '''
    multiprocessing.Array types are thread-safe by default, but are 
    horribly inefficient in getting/setting data. If you want speed you
    need to create a Numpy array pointing to the same shared memory,
    but this circumvents the automatic acquire/release behaviour. To 
    provide thread-safe behaviour would therefore require carrying through
    both the original Array (for the lock) and the derived Numpy array.
    This class is an attempt to get the best of both worlds: it behaves
    just like a normal Numpy array, but carries through the Array lock
    object and its methods. To use it:
    
    (In master process):
    import multiprocessing as mp
    mp_array = mp.Array(type, data_or_init)
    shared_numpy = SharedNumpyArray(mp_array)
    
      Pass shared_numpy to the thread Pool init function or to the thread
      itself if creating threads on the fly.
    
    (In each thread):
      If thread safety is not required (that is, different threads don't
      attempt to read and write to the same index), then just use it like
      any other array. If thread safety *is* required:
    
    with shared_numpy.get_lock():
        do_something(shared_numpy)
    '''
    
    def __new__(cls, mp_array):
        if mp_array is None:
            raise TypeError('Please provide a TypedMPArray object\
                             with a thread lock!')
        obj = frombuffer(mp_array.get_obj(), mp_array.dtype).view(cls)
        obj._mparray = mp_array
        obj.get_lock = mp_array.get_lock
        obj.acquire = mp_array.acquire
        obj.release = mp_array.release
        return obj

    def __array_finalize__(self, obj):
        if obj is None: 
            return
        
        self._mparray = getattr(obj, '_mparray', None)
        self.get_lock = getattr(obj, 'get_lock', None)
        self.acquire = getattr(obj, 'acquire', None)
        self.release = getattr(obj, 'release', None)
    
            
            
            
        
