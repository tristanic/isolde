import multiprocessing as mp
import ctypes
from . import SharedNumpyArray, TypedMPArray

class ChangeTracker:
    '''
    Handles communication between the master and slave threads via shared
    memory arrays.
    '''

    ALL_OUTPUTS           = 0xffff0000
    ALL_INPUTS            = 0x0000ffff

    def __init__(self):
        self._change_var = mp.Value(ctypes.c_int32, 0)
        self._managed_arrays = dict()
        self._possible_modified_arrays = dict()

    def clear_inputs(self):
        self._change_var.value &= ~self.ALL_INPUTS

    def clear_outputs(self):
        self._change_var.value &= ~self.ALL_OUTPUTS

    def run_all_necessary_callbacks(self, thread, changes):
        '''
        Runs the callbacks to apply the changes for any arrays whose
        contents have changed since the last iteration. Callbacks should
        be defined in the thread class, and take as arguments a Boolean array
        identifying which entries have changed, and a tuple of
        SharedNumpyArray objects containing the actual information.
        '''
        pma = self._possible_modified_arrays
        for change_bit in pma.keys():
            if changes & change_bit:
                for managed_array in pma[change_bit]:
                    key = managed_array[-2]
                    callback = getattr(thread, managed_array[-1])
                    change_mask, arrays = managed_array[2:4]
                    callback(change_mask, key, arrays)


    def add_managed_arrays(self, change_bit_mask, key, arrays, callback_name):
        '''
        Registers a combined set of shared arrays with the change
        tracker, and the name of the callback function that the
        simulation handler should run when it changes. Once registered,
        call register_array_changes(key, [change_mask]) to let the
        manager know that the array has changed, and optionally which
        elements have changed.

        @param change_bit_mask:
            The bit (e.g. ChangeTracker.DISTANCE_RESTRAINT) to be flipped
            to notify the top-level change tracker
        @param key:
            The key to the array in the SimComms dict
        @param arrays:
            A tuple containing the arrays themselves. All arrays must be
            the same length in the first dimension.
        @param callback_name:
            The name of the callback function in the SimThread to be
            run when the array changes.
        '''
        length = len(arrays[0])
        if len(arrays) > 1:
            for arr in arrays[2:]:
                if len(arr) != length:
                    raise TypeError(
                        'All arrays must be the same length in their '
                       +'first dimension!')

        change_flag = mp.Value(ctypes.c_bool, False)
        change_mask = SharedNumpyArray(TypedMPArray(ctypes.c_bool, length))


        ma = self._managed_arrays[key] =\
            (change_bit_mask, change_flag, change_mask, arrays, key, callback_name)
        try:
            self._possible_modified_arrays[change_bit_mask].append(ma)
        except KeyError:
            cl = self._possible_modified_arrays[change_bit_mask] = []
            cl.append(ma)

    def get_managed_arrays(self, key):
        '''
        Return the shared arrays associated with the given key.
        '''
        l = self._managed_arrays[key]
        return l[3]


    def register_array_changes(self, key, change_mask = None, indices = None):
        '''
        Let the change tracker know that an array has changed.
        @param key:
            The key to this array in the SimComms dict
        @param change_mask:
            A numpy Boolean array of the same length as the first
            dimension of the array, True where an element has changed.
            If neither change_mask nor indices is provided, it will be 
            assumed that all elements of the array have changed.
        @param indices:
            A numpy int array listing the indices that have changed.
            If neither change_mask nor indices is provided, it will be 
            assumed that all elements of the array have changed.
        '''
        if change_mask is not None and indices is not None:
            raise TypeError("Can't provide both change_mask and indices!")
        _change_bit, _change_flag, _change_mask = self._managed_arrays[key][0:3]
        with _change_flag.get_lock(), _change_mask.get_lock():
            _change_flag.value = True
            if change_mask is not None:
                numpy.logical_or(change_mask, _change_mask, out=_change_mask)
            elif indices is not None:
                _change_mask[indices] = True
            else:
                _change_mask[:] = True

        self.register_change(_change_bit)


    def clear_array_changes(self, key):
        _change_flag, _change_mask = self._managed_arrays[key][1:3]
        with _change_flag.get_lock(), _change_mask.get_lock():
            _change_flag.value = False
            _change_mask[:] = False

    @property
    def changes(self):
        return self._change_var

    def error(self):
        self.changes.value |= self.ERROR

    def register_change(self, change):
        changes = self.changes
        with changes.get_lock():
            changes.value |= change
