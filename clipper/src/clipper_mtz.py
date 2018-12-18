# @Author: Tristan Croll
# @Date:   28-Feb-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 29-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



#from . import clipper
from chimerax.core.models import Model

def calculate_voxel_size(resolution, shannon_rate):
    return resolution.limit/2/shannon_rate

def calculate_shannon_rate(resolution, voxel_size):
    return resolution.limit/(2*voxel_size)

class ReflectionDataContainer(Model):
    '''
    A container class to hold a set of reciprocal-space data, and
    defining the methods to access and use it. A sub-class of the
    ChimeraX Model class allowing it to be loaded into the model
    hierarchy making it easily visible to the user.
    '''
    def __init__(self, session, hklfile, shannon_rate = 1.5, min_voxel_size = 0.5,
        free_flag_label = None):
        '''
        This class should hold the information that's common to all
        the data contained in its children (e.g. the HKLinfo object,
        the Cell, Spacegroup and Grid_sampling objects, the Unit_Cell,
        etc.
        '''
        import os
        self.filename = os.path.basename(hklfile)
        Model.__init__(self, 'Reflection Data', session)
        hklinfo, free, exp, calc = load_hkl_data(session, hklfile, free_flag_label=free_flag_label)
        self._hklinfo = hklinfo
        self._grid_sampling = None

        voxel_size = calculate_voxel_size(hklinfo.resolution, shannon_rate)
        if voxel_size < min_voxel_size:
            shannon_rate = calculate_shannon_rate(hklinfo.resolution, min_voxel_size)
        self.shannon_rate = shannon_rate

        if free[0] is not None:
            self.free_flags = ReflectionData_FreeFlags(free[0], self.session, free[1])
            self.add([self.free_flags])

        self._experimental_data = None
        if len(exp[0]):
            dsets = []
            for name, data in zip(*exp):
                    dsets.append(ReflectionData_Exp(name, self.session, data))
            self.experimental_data.add(dsets)

        self._calculated_data = None
        if len(calc[0]):
            dsets = []
            for name, data in zip(*calc):
                    dsets.append(ReflectionData_Calc(name, self.session, data))
            self.calculated_data.add(dsets)

    @property
    def experimental_data(self):
        ed = self._experimental_data
        if ed is None or ed.deleted:
            ed = self._experimental_data=ReflectionData_Node('Experimental', self.session)
            self.add([ed])
        return ed

    @property
    def calculated_data(self):
        cd = self._calculated_data
        if cd is None or cd.deleted:
            cd = self._calculated_data=ReflectionData_Node('Calculated', self.session)
            self.add([cd])
        return cd

    @property
    def hklinfo(self):
        return self._hklinfo

    @property
    def cell(self):
        return self.hklinfo.cell

    @property
    def spacegroup(self):
        return self.hklinfo.spacegroup

    @property
    def resolution(self):
        return self.hklinfo.resolution

    @property
    def grid_sampling(self):
        if self._grid_sampling is None:
            from . import Grid_sampling
            self._grid_sampling = Grid_sampling(
                self.spacegroup, self.cell, self.resolution, self.shannon_rate)
        return self._grid_sampling





class ReflectionData_Node(Model):
    '''
    Container class to hold a subset of reflection data within a
    ReflectionDataContainer tree. Typically the subset will be either
    'Experimental' or 'Calculated'.
    '''
    # def __init__(self, name, session):
    #     Model.__init__(self, name, session)
    #     # self.datasets = datasets
    #     for name, data in datasets.items():
    #         self.add([data])

    @property
    def datasets(self):
        return dict((m.name, m) for m in self.child_models())

    def __iter__(self):
        return iter(self.datasets.values())

    def __getitem__(self, key):
        return self.datasets[key]

    def __len__(self):
        return len(self.datasets)

class ReflectionData(Model):
    '''
    Prototype for ReflectionData_Exp and ReflectionData_Calc. Should
    contain methods common to both (e.g. drawing of reciprocal-space
    reflections).
    '''
    def __init__(self, name, session, data):
        '''
        Args:
            name:
                A descriptive name.
            session:
                The ChimeraX session.
            data:
                A Clipper HKL_data_Flag or HKL_data_Flag_bool object.
        '''

        Model.__init__(self, name, session)
        self._data = data

    @property
    def data(self):
        return self._data

    @property
    def dtype(self):
        return type(self.data)



class ReflectionData_FreeFlags(ReflectionData):
    '''Holds the array of free flags.'''
    pass


class ReflectionData_Exp(ReflectionData):
    '''
    Holds one set of experimental reflection data (e.g. F/sigF, I/sigI),
    etc.
    '''
    pass

class ReflectionData_Calc(ReflectionData):
    '''Holds one set of calculated reflections and phases.'''
    def __init__(self, name, session, data, is_difference_map = None):
        '''
        Args:
            name:
                A descriptive name.
            session:
                The ChimeraX session.
            data:
                A Clipper HKL_data_F_phi object
            is_difference_map(bool):
                If True, maps generated from this data will be displayed
                with both positive and negative contours. If None, the
                code will attempt to guess from the name
        '''
        ReflectionData.__init__(self, name, session, data)
        if is_difference_map is None:
            self.is_difference_map = self._guess_if_difference_map(name)
        else:
            self.is_difference_map = is_difference_map

    def _guess_if_difference_map(self, name):
        if name.split(',')[0] in ('F','FWT'):
            return False
        if '2' in name:
            return False
        return True


def data_labels_match(data_col, sigma_or_phase_col):
  '''
  A quick-and-dirty approach to check if MTZ column names belong together.
  Relies on the fact that under most standard naming the phi or sigma column
  name corresponds to the amplitude/intensity column plus a prefix. We'll
  pair together columns where the second one matches the first plus a prefix
  or a suffix, but not both. We also need to make sure that we don't
  accidentally pair (for example) FOFC with PH2FOFC, so we'll exclude
  any cases where the prefix/suffix contains a digit
  '''
  if data_col in sigma_or_phase_col:
    remainder = [s for s in sigma_or_phase_col.split(data_col) if s]
    if len(remainder) == 1:
      return not any(i.isdigit() for i in remainder[0])
  return False

def find_data_pairs(mtzin, temp_tree, first_type, second_type, data_type):
    '''
    Find all column pairs matching a particular signature and with
    matching labels, create the corresponding Clipper objects, and
    return the objects and their names as two arrays.

    Args:
        mtzin:
            A currently open clipper.CCP4MTZfile object
        temp_tree:
            A DataTree (or similar dict) mirroring the internal data
            structure of the MTZ file.
        first_type (char):
            a character defining the data type for the first column. See:
            http://www.ccp4.ac.uk/html/mtzformat.html#coltypes
            for a list of valid types.
        second_type (char):
            a character defining the data type for the second column.
        data_type:
            a pointer to the clipper class corresponding to this
            data type.
    '''
    data_sets = []
    data_set_names = []
    for crystal in temp_tree.keys():
        for dataset in temp_tree[crystal].keys():
            for iname in temp_tree[crystal][dataset][first_type].keys():
                for sname in temp_tree[crystal][dataset][second_type].keys():
                    if data_labels_match(iname, sname):
                        keyname = ', '.join(map(str, [iname, sname]))
                        this_set = data_type()
                        mtzin.import_hkl_data(this_set,
                            '/'+'/'.join(
                            map(str, [crystal, dataset, '[{}]'.format(keyname)])))
                        data_sets.append(this_set)
                        data_set_names.append(keyname)
    return (data_sets, data_set_names)

def find_free_set(mtzin, temp_tree, label = None):
    '''Find the free set, optionally given an explicit label to look for.'''
    possible_free_flags = []
    possible_free_flags_names = []
    from . import HKL_data_Flag
    for crystal in temp_tree.keys():
        for dataset in temp_tree[crystal].keys():
            # Find and store the free set
            for iname in temp_tree[crystal][dataset]['I'].keys():
                if label is not None:
                    if iname == label:
                        free = HKL_data_Flag()
                        mtzin.import_hkl_data(free,
                            '/'+'/'.join(map(str, [crystal, dataset, '[{}]'.format(iname)])))
                        return ([free],[iname])

                elif 'free' in iname.lower():
                    thisFree = HKL_data_Flag()
                    mtzin.import_hkl_data(thisFree,
                        '/'+'/'.join(map(str, [crystal, dataset, '[{}]'.format(iname)])))
                    possible_free_flags.append(thisFree)
                    possible_free_flags_names.append(iname)

    if label is not None:
        raise TypeError('The label "{}" does not appear to be in this dataset!'.format(label))
    return (possible_free_flags, possible_free_flags_names)


def load_hkl_data(session, filename, free_flag_label = None):
    '''
    Load in an mtz file, create Clipper objects from the data, and
    return the tuple:
        ( HKLinfo,
         (free_flags_name, free_flags),
         (experimental_set_names, experimental_sets),
         (calculated_set_names, calculated_sets) )

    where HKLinfo and free_flags are Clipper objects, and
    experimental_sets and calculated_sets are arrays of Clipper objects.

    If free_flag_label is a string, the column with that name will be assigned
    as the free flags. If free_flag_label is -1, we will continue even if
    no free set is found. If multiple possible free sets are found, or
    no free sets are found and free_flag_label != -1, an error will be
    raised.
    '''
    import os
    if os.path.isfile(filename):
        hklfile = os.path.abspath(filename)
    else:
        raise FileNotFoundError('Invalid filename!')

    extension = os.path.splitext(filename)[1].lower()
    if extension not in ('.cif', '.mtz'):
        raise ValueError('Reflection file must be in either .cif or .mtz format!')
        return

    if extension == '.mtz':
        from . import CCP4MTZfile, HKL_info
        mtzin = CCP4MTZfile()
        hkl = HKL_info()
        mtzin.open_read(filename)
        mtzin.import_hkl_info(hkl)
        # Get all the column names and types
        column_labels = mtzin.column_paths
        # Sort the columns into groups, and organise into a temporary tree
        from .data_tree import DataTree
        temp_tree = DataTree()['Experiment']
        i = 0
        for l in column_labels:
            thisname, thistype = l.__str__().split(' ')
            crystal, dataset, name = thisname.split('/')[1:]
            # The h, k, l indices are already captured in self.hklinfo
            if thistype != 'H':
                temp_tree[crystal][dataset][thistype][name] = i
                i += 1

        # FIXME: will need checking here to see if the crystal name already
        # exists, and offer options to replace, rename or append

        # FIXME: ultimately we want to be able to handle any legal MTZ file,
        #        but for now we'll balk if we find more than one crystal
        if len(temp_tree.keys()) != 1:
            errstring =\
            '''
            At present ChimeraX-Clipper cannot handle MTZ files containing
            data from more than one crystal. We hope to fix this in a future
            release, but for now you can split your MTZ file using tools
            available in the Phenix or CCP4 suites. Note that you *can*
            load multiple MTZ files into the one Clipper_MTZ datastructure
            if you wish, using repeat calls to load_hkl_data. However, they
            will all be treated as having the same unit cell parameters. If
            in doubt, create a new Xtal_Project object for each file.
            '''
            raise RuntimeError(errstring)

        # Find experimental data sets
        experimental_sets = []
        experimental_set_names = []


        # Find I/sigI pairs and pull them in as
        # Clipper HKL_data_I_sigI objects
        from . import HKL_data_I_sigI
        from . import HKL_data_F_sigF
        from . import HKL_data_F_phi
        isigi, isigi_names = find_data_pairs(
            mtzin, temp_tree, 'J', 'Q', HKL_data_I_sigI)
        experimental_sets.extend(isigi)
        experimental_set_names.extend(isigi_names)


        # Find F/sigF pairs and pull them in as
        # Clipper HKL_data_F_sigF objects
        fsigf, fsigf_names = find_data_pairs(
            mtzin, temp_tree, 'F', 'Q', HKL_data_F_sigF)
        experimental_sets.extend(fsigf)
        experimental_set_names.extend(fsigf_names)

        calculated_sets = []
        calculated_set_names = []

        # Find amplitude/phase pairs and pull them in as
        # Clipper HKL_data_F_phi objects
        fphi, fphi_names = find_data_pairs(
            mtzin, temp_tree, 'F', 'P', HKL_data_F_phi)
        calculated_sets.extend(fphi)
        calculated_set_names.extend(fphi_names)

        free_flags = None
        free_flags_name = None

        possible_free_flag_names = []
        all_integer_column_names = []
        for crystal in temp_tree.keys():
            for dataset in temp_tree[crystal].keys():
                for iname in temp_tree[crystal][dataset]['I'].keys():
                    if 'free' in iname.lower():
                        possible_free_flag_names.append(iname)
                    all_integer_column_names.append(iname)

        # possible_free_flags, possible_free_flag_names = find_free_set(
        #     mtzin, temp_tree, free_flag_label)

        if (len(possible_free_flag_names) == 0 and len(experimental_sets)):
            all_integer_column_names = temp_tree[crystal][dataset]['I'].keys()
            if not len(all_integer_column_names) and free_flag_label !=-1:
                err_string = 'This MTZ file does not appear to contain any '\
                + 'columns suitable for use as a free set. Please generate '\
                + 'a suitable set of free flags using your favourite '\
                + 'package (I recommend PHENIX) and try again.'
                raise RuntimeError(err_string)
            else:
                possible_free_flag_names = all_integer_column_names
                if len(possible_free_flag_names) == 1:
                    session.logger.info('WARNING: assuming column with label {}'
                    + ' defines the free set.'.format(possible_free_flag_names[0]))
                    #possible_freefrom ._flags = [temp_tree[crystal][dataset]['I'][possible_free_flag_names[0]]]
        if not possible_free_flag_names:
            free_flags_name = None
        elif len(possible_free_flag_names) > 1:
            free_flags_name = _r_free_chooser(session, possible_free_flag_names)
            if free_flags_name is None:
                if free_flag_label != -1:
                    raise RuntimeError('No free flags chosen. Bailing out.')
        else:
            free_flags_name = possible_free_flag_names[0]

        if free_flags_name:
            free_flags = find_free_set(mtzin, temp_tree, label=free_flags_name)[0][0]
        else:
            free_flags = None

        mtzin.close_read()

    return ( hkl,
            (free_flags_name, free_flags),
            (experimental_set_names, experimental_sets),
            (calculated_set_names, calculated_sets) )

def _r_free_chooser(session, possible_names):
    from PyQt5.QtWidgets import QInputDialog
    choice, ok_pressed = QInputDialog.getItem(session.ui.main_window, 'Choose R-free column', 'Label: ', possible_names, 0, False)
    if ok_pressed and choice:
        return choice
    return None
