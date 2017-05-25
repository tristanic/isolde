import numpy
from chimerax.core.atomic import Residue, AtomicStructure
def is_continuous_protein_chain(sel):
    '''
    Checks if the residues in a selection are all protein, and form a 
    continuous chain with no breaks. 
    NOTE: only one atom from each residue need be selected.
    '''
    try:
        p = find_polymer(sel)
    except:
        return False
    res = sel.unique_residues
    if not numpy.all(res.polymer_types == Residue.PT_AMINO):
        return False
        
    indices = sorted(p.indices(res))
    first = indices[0]
    # Check if the indices are continuous
    return all(a == b for a, b in enumerate(indices, indices[0]))    

def find_polymer(sel):
    '''
    Find the polymer containing the first residue in a selection.
    '''
    us = sel.unique_structures
    if len(us) != 1:
        raise TypeError('All atoms must be from the same structure!')
    m = us[0]
    res = sel.unique_residues
    polymers = m.polymers(
        missing_structure_treatment = m.PMS_NEVER_CONNECTS)
    first_index = -1
    r0 = res[0]
    for p in polymers:
        first_index = p.index(r0)
        if first_index != -1:
            return p
    raise IndexError('Polymer not found!')
            
   
def add_disulfides_from_model_metadata(model):
    metadata = model.metadata
    
