# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tcroll@altoslabs.com
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

def expand_selection_to_neighbors(atoms, pad):
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
    sel = sm.atoms[sm.atoms.selected]
    sel.intra_bonds.selected = True
    return sel


def extend_selection_along_chains(residues, direction):
    '''
    Extends a selection by one residue in the given direction, stopping
    when it hits a chain break or the end of a chain. Every contiguous
    stretch of polymeric residues will be extended in this fashion. It 
    is expected that all residues are from the same structure.
    Args:
        residues:
            A ChimeraX :class:`Residues` object
        direction:
            -1 or 1
    '''
    from chimerax.atomic import Residues, concatenate
    new_residues = []
    m = residues.unique_structures[0]
    polymers = m.polymers(
        missing_structure_treatment = m.PMS_NEVER_CONNECTS)
    def _extend(p, frag):
        if direction==-1:
            prev_index = frag[0]-1
            if prev_index >= 0:
                new_residues.append(p[prev_index])
        else:
            last_index = frag[-1]
            next_index = last_index+1
            if next_index < len(p):
                new_residues.append(p[next_index])
    for p in polymers:
        p = p[0]
        indices = p.indices(residues)
        indices = indices[indices!=-1]
        if not len(indices):
            continue
        indices = list(sorted(indices))
        while len(indices):
            last_index = indices[0]
            frag = [last_index]
            for i, index in enumerate(indices[1:]):
                if index-last_index==1:
                    frag.append(index)
                    last_index = index
                else:
                    _extend(p, frag)
                    indices = indices[i+1:]
                    break
            else:
                _extend(p, frag)
                break 
    residues = concatenate([residues, Residues(new_residues)])
    residues.atoms.selected=True
    residues.atoms.intra_bonds.selected=True

    return residues


def shrink_selection_by_one_res(residues, direction):
    '''
    For each set of contiguous polymeric residues, shrinks the current selection 
    by one residue from one end. If direction == -1 the first residue will be removed 
    from the selection, otherwise the last will be removed. Each contiguous selection 
    will never be shrunk to less than a single residue. It is expected that all residues 
    are from the same structure.
    '''
    
    m = residues.unique_structures[0]
    for p in m.polymers(missing_structure_treatment=m.PMS_NEVER_CONNECTS):
        p = p[0]
        indices = p.indices(residues)
        indices = indices[indices!=-1]
        if not len(indices):
            continue
        indices = list(sorted(indices))
        def _shrink(p, frag):
            if len(frag)==1:
                return
            if direction==-1:
                r = p[frag[0]]
            else:
                r = p[frag[-1]]
            r.atoms.selected=False
            r.atoms.intra_bonds.selected=False
            for n in r.neighbors:
                r.bonds_between(n).selected=False
        while len(indices):
            last_index = indices[0]
            frag = [last_index]
            for i, index in enumerate(indices[1:]):
                if index-last_index==1:
                    frag.append(index)
                    last_index = index
                else:
                    _shrink(p, frag)
                    indices = indices[i+1:]
                    break
            else:
                _shrink(p, frag)
                break