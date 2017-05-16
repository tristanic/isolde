import numpy
from chimerax.core.atomic import Residue
def is_continuous_protein_chain(sel):
    '''
    Checks if the residues in a selection are all protein, and form a 
    continuous chain with no breaks. 
    NOTE: only one atom from each residue need be selected.
    '''
    us = sel.unique_structures
    if len(us) != 1:
        return False
    m = us[0]
    res = sel.unique_residues
    if not numpy.all(res.polymer_types == Residue.PT_AMINO):
        return False
    polymers = m.polymers(consider_missing_structure = False)
    first_index = -1
    r0 =res[0]
    for p in polymers:
        first_index = p.index(r0)
        if first_index != -1:
            break
    
    if first_index == -1:
        return False
    
    indices = sorted(p.indices(res))
    first = indices[0]
    # Check if the indices are continuous
    return all(a == b for a, b in enumerate(indices, indices[0]))    
            
    
    
