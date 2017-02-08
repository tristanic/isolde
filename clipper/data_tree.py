'''
Main organisational tree of crystal data. A typical MTZ file
is ordered like:

    /HKL_base/HKL_base/H H
    /HKL_base/HKL_base/K H
    /HKL_base/HKL_base/L H
    /crystal/dataset/2FOFCWT F
    /crystal/dataset/PH2FOFCWT P
    /crystal/dataset/2FOFCWT_fill F
    /crystal/dataset/PH2FOFCWT_fill P
    /crystal/dataset/FOFCWT F
    /crystal/dataset/PHFOFCWT P

... where the final digit defines the data type, and the leading
path allows the file to hold data from multiple experiments/crystals.
The HKL coordinate arrays (type 'H') are first split off into a 
HKL_info object, and the remaining data is ordered first by the two 
levels of the path, followed by type and finally by name. 
'''

import collections
from enum import IntEnum
class db_levels(IntEnum):
  PROJECT = 0
  EXPERIMENT = 1
  CRYSTAL = 2
  DATASET = 3


class DataTree(collections.defaultdict):
  '''
  A simple extension to the defaultdict object, which automatically
  builds a tree data structure. Contains a dict to hold metadata.
  Metadata for a given node should apply to that node and all its
  children. When metadata is requested from a node, it searches for
  the key first in its own _metadata dict, and if not found passes the
  request to its parent, and so on. A typical data tree (and the 
  metadata at each level) might be:
  
  [Project] (name, Lab head, Experimenter, Detailed description, ...)
     |
     -- [Experiment] (Name, protein(s), ligand(s), ...)
     |      |
     |      -- [Crystal] (Name, conditions, unit cell, ...)
     .      |     |
     .      |     -- [Dataset] (Name, wavelength, beamline, ...)
     .      .     |   |
            .     |   -- HKL_info
            .     .   |
                  .   -- HKL_data_F_sigF
                  .   |
                      -- HKL_data_F_sigF_ano
                      |
                      -- HKL_data_F_phi
                      |
                      -- Xmap
                      |
                      -- etc.
  ''' 
     
  def __init__(self, *args, parent = None, project = 'Project', experiment = 'Experiment'):
    collections.defaultdict.__init__(self, *args)
    # Keep track of our parents...
    self.parent = parent
    # ... and what level we are at in the tree
    if self.parent is None:
        self._level = db_levels.PROJECT
        self._metadata = {'Name': project}
        # Initialise the first experiment node
        self[experiment].set_metadata(('Name', experiment))
    else:
        self._level = self.parent.level + 1
    self._metadata = {}
  
  @property
  def level(self):
      return self._level
      
  def find_ancestor(self, level = db_levels.PROJECT):
    '''
    Return the node at the given level within the tree.
    '''
    if self.level == level:
      return self
    return self.parent.find_ancestor(level)
  
  def set_metadata(self, *key_data_pairs):
    for key, arg in key_data_pairs:
      self._metadata[key] = arg
  
  def get_metadata(self, key):
    try:
      return self._metadata[key]
    except KeyError:
      if self.parent is not None:
        return self.parent.get_metadata(key)
      else:
        raise
  
  def __missing__(self, key):
    ret = self[key] = DataTree(parent = self)
    ret.set_metadata(('Name', key))
    return ret
