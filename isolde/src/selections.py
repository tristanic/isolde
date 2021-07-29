# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



def get_shell_of_residues(residues, dist_cutoff):
    '''
    Get a shell of whole residues from the same model as the atoms in residues,
    within a user-defined cut-off distance surrounding residues. Expects
    all residues to be within the same model.
    '''
    from chimerax.geometry import find_close_points
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

def expand_selection_along_chains(atoms, pad):
    '''
    Expand an existing selection to whole residues, and extend outwards by `pad` covalently-bonded 
    neighbors.
    '''
    us = atoms.unique_structures
    if len(us) != 1:
        raise TypeError('Selected atoms must all be in the same model!')
    sm = us[0]
    session = sm.session
    selatoms = sm.atoms[sm.atoms.selected]
    from chimerax.atomic import selected_atoms
    sm.atoms.selected=False
    sm.bonds.selected=False
    selres = selatoms.unique_residues
    from .util import expand_selection
    expand_selection(selres, pad)
    return sm.atoms[sm.atoms.selected]
