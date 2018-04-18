# @Author: Tristan Croll
# @Date:   22-Mar-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   Tristan Croll
# @Last modified time: 18-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



def get_shell_of_residues(residues, dist_cutoff):
    '''
    Get a shell of whole residues from the same model as the atoms in residues,
    within a user-defined cut-off distance surrounding residues. Expects
    all residues to be within the same model.
    '''
    from chimerax.core.geometry import find_close_points
    from chimerax.atomic import selected_atoms, Atoms, concatenate
    us = residues.unique_structures
    selatoms = residues.atoms
    if len(us) !=1:
        raise Exception('selection should contain atoms from a single molecule!')
    allres = us[0].residues
    unsel_residues = allres.subtract(residues)
    unselected_atoms = unsel_residues.atoms
    selcoords = selatoms.coords
    unselcoords = unselected_atoms.coords
    ignore, shell_indices = find_close_points(selcoords, unselcoords, dist_cutoff)
    shell_residues = unselected_atoms[shell_indices].unique_residues
    return shell_residues

def expand_selection_along_chains(atoms, extension):
    '''
    Expand an existing selection to whole residues, and extend backwards and
    forwards along each chain by the number of residues defined by extension,
    stopping at chain breaks.
    '''
    us = atoms.unique_structures
    if len(us) != 1:
        raise TypeError('Selected atoms must all be in the same model!')
    m = us[0]
    residues = atoms.unique_residues
    all_residues = m.residues
    res_indices = all_residues.indices(residues)

    from chimerax.atomic import Structure
    all_frags = m.polymers(missing_structure_treatment = Structure.PMS_NEVER_CONNECTS)
    import numpy
    sel_frags = []
    sel_frag_res_indices = []
    master_sel_mask = numpy.zeros(len(all_residues), numpy.bool)
    for frag in all_frags:
        frag = frag[0]
        if numpy.any(frag.indices(residues) > -1):
            sel_frags.append(frag)
            sel_frag_res_indices.append(all_residues.indices(frag))
    for frag, frag_indices in zip(sel_frags, sel_frag_res_indices):
        frag_nres = len(frag_indices)
        sel_mask = numpy.isin(frag_indices, res_indices, assume_unique=True)
        sel_pol_indices = numpy.where(sel_mask)[0]
        for i in sel_pol_indices:
            lb = i-extension
            ub = i+extension+1
            if lb<0:
                lb = 0
            if ub > frag_nres:
                ub = frag_nres
            sel_mask[lb:ub] = True
        master_sel_mask[frag_indices[sel_mask]] = True
    return all_residues[master_sel_mask].atoms
