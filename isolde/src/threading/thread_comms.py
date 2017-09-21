import multiprocessing as mp

class ThreadComms(dict):
    '''
    Holds all the shared variables/arrays needed to communicate with a
    thread.
    '''

    def __init__(self):
        self._locked = False
        super().__init__()

    def __setitem__(self, key, obj):
        if self._locked:
            raise KeyError('Cannot change the composition of a {} object '\
                +'while a thread is running!'.format(self.__class__))
        super().__setitem__(key, obj)

    def lock(self):
        '''
        Lock the comms object to prevent addition/removal of entries while a
        thread is running.
        '''
        self._locked = True

    def unlock(self):
        '''
        Unlock the comms object to allow all modifications. This should only
        be done while it is not in use by any child threads.
        '''
        self._locked = False

    def thread_safe_set_value(self, key, val):
        '''Set the value of a mp.Value object in a thread-safe way.'''
        target = self[key]
        with target.get_lock():
            target.value = val

    def thread_safe_set_array_values(self, key, values, indices_or_mask = None):
        '''Change values within a SharedNumpyArray in a thread-safe way.'''
        target = self[key]
        with target.get_lock():
            if indices_or_mask is not None:
                target[indices_or_mask] = values
            else:
                target[:] = values
