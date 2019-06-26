# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 21-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll



def _get_metals():
    '''
    Returns a list of metals sorted by mass
    '''
    from chimerax.atomic import Element
    metals = [n for n in Element.names if Element.get_element(n).is_metal]
    sorted_metals = sorted(metals, key=lambda m: Element.get_element(m).mass)
    return [Element.get_element(m) for m in sorted_metals]

_metals = _get_metals()

def place_metal_at_coord(model, chain_id, residue_number, residue_name, atom_name, coord, element_name=None, bfactor=20):
    '''
    Create a new residue encompassing a single metal ion in the current model,
    and place it at the given coordinate.
    '''
    if element_name is None:
        element_name = atom_name.title()
    from chimerax.atomic import Element
    e = Element.get_element(element_name)
    r = model.new_residue(residue_name, chain_id, residue_number)
    a = model.new_atom(atom_name, e)
    a.coord = coord
    a.bfactor=bfactor
    r.add_atom(a)

def find_nearest_chain(model, coord):
    '''
    Find the chain coming nearest to the given coordinate.
    '''
    pass

def suggest_new_residue_number_for_ligand(chain):
    '''
    Suggest a suitable residue number for a new ligand, based on what is already
    in the chain.
    '''
    pass
