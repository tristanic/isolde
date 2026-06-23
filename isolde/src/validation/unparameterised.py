# @Author: Tristan Croll <tic20>
# @Date:   01-Aug-2020
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 28-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

# Make sure core CIF templates are pre-loaded
from chimerax.isolde import atomic

H_TO_HEAVY_ATOM_THRESHOLD_RATIO = 0.5

def suspiciously_low_h(residues):
    '''
    Heuristic: returns True if the H/heavy-atom ratio of the given residues is
    below ``H_TO_HEAVY_ATOM_THRESHOLD_RATIO``. Always False if there are no
    heavy atoms.
    '''
    enames = residues.atoms.element_names
    n_h = (enames == 'H').sum()
    n_heavy = (enames != 'H').sum()
    if n_heavy == 0:
        return False
    return (n_h / n_heavy) < H_TO_HEAVY_ATOM_THRESHOLD_RATIO

def waters_without_h(residues):
    '''
    Returns True if any residue named 'HOH' in the given set has other than
    three atoms (i.e. is missing one or both hydrogens).
    '''
    waters = residues[residues.names == 'HOH']
    for w in waters:
        if len(w.atoms) != 3:
            return True
    return False
