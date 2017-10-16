import numpy

class StructureStepper:
  '''
  Class to aid in stepping through an atomic model some number of
  secondary structure elements (plus flanking loops) at a time.
  '''
  def __init__(self, session, model, step = 2, min_gap = 2):
    '''
    Args:
      session:
        The ChimeraX session
      model:
        An AtomicStructure object
      step (int):
        The number of secondary structure elements to be covered in each
        step
      min_gap (int):
        If a secondary structure element contains fewer than this number
        of residues, it will be rolled into the next element
    '''
    self.session = session
    self.model = model  
    self.step = step
    self.min_gap = min_gap
    self._num_polymers = len(model.polymers()) 
    cp = self.current_polymer = model.polymers()[0]
    self.min_element = cp[0].secondary_structure_ids[0]
      
  def step_forward(self):
    '''
    Step through the structure a given number of defined secondary structure
    elements at a time. Not yet perfect, but working. Needs to be adapted
    to allow stepping forward and backward, extending the selection, etc.
    Also needs methods for stepping through other polymers (particularly
    nucleic acids).
    '''
    from chimerax.core.atomic import Residue
    min_element = self.min_element
    cres, ptype = self.current_polymer
    ptype = self.current_polymer[1]
    done = False
    chain_finished = False
    all_finished = False
    element_count = 0
    this_element = min_element
    last_type = None
    if cres[0].polymer_type == Residue.PT_AMINO:
      # Step through by secondary structure
      res_indices = numpy.where(cres.secondary_structure_ids == this_element)[0]
      residues = cres[res_indices]
      while not done:
        if not len(residues):
          #We've passed the end of the chain. Next time we start on the
          # next chain (if it exists)
          chain_finished = True
          if self.model.polymers().index(self.current_polymer) == self._num_polymers - 1:
            all_finished = True
          break
        r = residues[0]
        this_type = _res_ss_type(residues[0])
        print('This element: {}, Current type: {}'.format(this_element, this_type))
        if len(residues) < self.min_gap:
          next_index = max(res_indices) + 1
          if next_index >= len(cres):
            chain_finished = True
            # End of chain.
            break
          next_type = _res_ss_type(cres[next_index])
          if next_type == last_type:
            # We assume the next element is an extension of the last
            element_count = max (element_count -1, 0)
          elif element_count >= self.step:
            # We're done
            break
        elif this_type == Residue.SS_STRAND or this_type == Residue.SS_HELIX:
          last_found_element = this_element
          element_count += 1
        print(element_count)
        this_element += 1
        last_type = this_type
        res_indices = numpy.where(cres.secondary_structure_ids == this_element)[0]
        residues = cres[res_indices]
        if len(residues) > self.min_gap and element_count >= self.step:
          done = True
      self.min_element = last_found_element
      ret = cres[numpy.logical_and(cres.secondary_structure_ids >= min_element, 
                          cres.secondary_structure_ids <= this_element)].atoms
      if chain_finished:
        if all_finished:
          self.session.logger.info("Reached the end of the model. Returning to start...")
          self.current_polymer = self.model.polymers()[0]
        else:
          pol = self.current_polymer = self.model.polymers()[self.model.polymers().index(self.current_polymer) + 1]
        self.min_element = pol[0].secondary_structure_ids[0]
      return ret
    else:
      # Simply step through in blocks of n residues.
      pass
    
      
      
def _res_ss_type(residue):
  from chimerax.core.atomic import Residue
  if residue.is_strand:
    return Residue.SS_STRAND
  if residue.is_helix:
    return Residue.SS_HELIX
  return -1
