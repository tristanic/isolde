from .data_tree import DataTree
from . import clipper
class Clipper_MTZ(DataTree):
  '''
  A container class to capture all the information loaded in from a 
  .mtz or .cif structure factor file.
  structure, using the clipper-python libraries: 
    - The reflection data
    - The cell and symmetry parameters
    - The map(s)
    - The atomic model describing the asymmetric unit
  
  Contains methods for:
    - Loading in reflection and reciprocal-space map data from 
      structure factor (.mtz or .cif) files
    - Displaying density surrounding the current center of rotation
    - Quickly finding and displaying all local symmetry equivalents
    ...
    
  NOTE: The modern MTZ file has the capability to hold data for multiple
  crystals of the same spacegroup, but with different cell dimensions.
  Officially, the global CELL definition is deprecated in favour of the
  crystal-specific DCELL, but Clipper still appears to only have functions
  to load in the single global definition. Hence, the below code will
  fail on such multi-crystal MTZ files. Luckily, in practice very few
  people use MTZ files in this way, and tools available in the Phenix
  or CCP4 suites can be used to split them into individual datasets. For
  now, we'll just raise an exception if we find more than one crystal in
  the file.
  '''

  def __init__(self, parent = None):
    '''
    parent can be either None or an Experiment-level node of an existing
    data tree. If None, an empty tree will be created to act as the parent
    '''
    if parent is None:
      parent = DataTree()['Experiment']
    DataTree.__init__(self, parent = parent)
    self.parent = parent
    self._temp_tree = None

  def has_map_data(self, crystal_name):
    if crystal_name in self.keys():
      for key, dataset in self[crystal_name].items():
        for key, data in dataset.items():
          if key == 'F_Phi' and len(data):
            return True
          else:
            return False
    else:
      raise KeyError('No such crystal!')
    
  def load_hkl_data(self, filename):
    '''
    Load in an mtz file, create Clipper objects from the data, and add
    them to the database. Returns the name of the crystal found in the
    mtz file.
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
    
    from datetime import datetime
    self.set_metadata(('Source file', hklfile))
    self.set_metadata(('Accessed',str(datetime.today())))
    
    if extension == '.mtz':
      mtzin = clipper.CCP4MTZfile()
      hkl = self['hkl'] = clipper.HKL_info()
      mtzin.open_read(filename)
      mtzin.import_hkl_info(hkl)
      # Get all the column names and types
      column_labels = mtzin.column_labels
      # Sort the columns into groups, and organise into a temporary tree 
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
        print(errstring)
        print(temp_tree)
        self._temp_tree = temp_tree
        return
      
      for crystal in temp_tree.keys():
        for dataset in temp_tree[crystal].keys():
          # Find and store the free set
          for iname in temp_tree[crystal][dataset]['I'].keys():
            if 'free' in iname.lower():
              assert 'Free flags' not in self[crystal][dataset].keys()
              thisFree = self[crystal][dataset]['Free flags'] \
                   = clipper.HKL_data_Flag()
              mtzin.import_hkl_data(thisFree,
                '/'+'/'.join(map(str, [crystal, dataset, '[{}]'.format(iname)])))
                
          # Find I/sigI pairs and pull them in as
          # Clipper HKL_data_I_sigI objects
          for iname in temp_tree[crystal][dataset]['J'].keys():
            for sname in temp_tree[crystal][dataset]['Q'].keys():
              if data_labels_match(iname, sname):
                keyname = ', '.join(map(str, [iname, sname]))
                thisIsigI \
                  = self[crystal][dataset]['Intensities'][keyname] \
                  = clipper.HKL_data_I_sigI()
                mtzin.import_hkl_data(thisIsigI,
                '/'+'/'.join(map(str, [crystal, dataset, '[{}]'.format(keyname)])))
          # Find F/sigF pairs and pull them in as
          # Clipper HKL_data_F_sigF objects
          for fname in temp_tree[crystal][dataset]['F'].keys():
            for sname in temp_tree[crystal][dataset]['Q'].keys():
              if data_labels_match(fname, sname):
                keyname = ', '.join(map(str, [fname, sname]))
                thisFsigF \
                  = self[crystal][dataset]['F_SigF'][keyname] \
                  = clipper.HKL_data_F_sigF()
                mtzin.import_hkl_data(thisFsigF, 
                '/'+'/'.join(map(str, [crystal, dataset, '[{}]'.format(keyname)])))
          
          # Find amplitude/phase pairs and pull them in as 
          # Clipper HKL_data_F_phi objects
          for fname in temp_tree[crystal][dataset]['F'].keys():
            for pname in temp_tree[crystal][dataset]['P'].keys():
              if data_labels_match(fname, pname):
                keyname = ', '.join(map(str, [fname, pname]))
                thisFphi \
                  = self[crystal][dataset]['F_Phi'][keyname] \
                  = clipper.HKL_data_F_phi()
                mtzin.import_hkl_data(thisFphi, 
                '/'+'/'.join(map(str, [crystal, dataset, '[{}]'.format(keyname)])))
    return crystal

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
