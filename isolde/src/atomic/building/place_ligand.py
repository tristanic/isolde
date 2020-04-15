# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 15-Apr-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

from . import set_new_atom_style


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

def place_water(session, model, position):
    from chimerax.core.geometry import find_closest_points
    from chimerax.atomic import mmcif
    matoms = model.atoms
    _,_,i = find_closest_points([position], matoms.coords, 3)
    na = matoms[i[0]]
    cid = na.residue.chain_id
    bfactor = na.residue.atoms.bfactors.mean()+5
    hoh = mmcif.find_template_residue(session, 'HOH')
    r = new_residue_from_template(model, hoh, cid, position, b_factor=    bfactor)
    matoms.selected=False
    r.atoms.selected=True
    from chimerax.core.commands import run
    run(session, 'isolde sim start sel')

def new_residue_from_template(model, template, chain_id, center,
        residue_number=None, insert_code=' ', b_factor=50, precedes=None):
    '''
    Create a new residue based on a template, and add it to the model.
    '''
    if residue_number is None:
        if chain_id in model.residues.chain_ids:
            residue_number = suggest_new_residue_number_for_ligand(model, chain_id)
        else:
            residue_number = 0
    import numpy
    from chimerax.atomic import Atom
    t_coords = numpy.array([a.coord for a in template.atoms])
    t_center = t_coords.mean(axis=0)
    t_coords += numpy.array(center) - t_center
    tatom_to_atom = {}
    r = model.new_residue(template.name, chain_id, residue_number,
        insert=insert_code, precedes=precedes)
    for i, ta in enumerate(template.atoms):
        a = tatom_to_atom[ta] = model.new_atom(ta.name, ta.element)
        a.coord = t_coords[i]
        a.bfactor = b_factor
        r.add_atom(a)
        for tn in ta.neighbors:
            n = tatom_to_atom.get(tn, None)
            if n is not None:
                model.new_bond(a, n)
    set_new_atom_style(model.session, r.atoms)
    return r



def find_nearest_chain(model, coord):
    '''
    Find the chain coming nearest to the given coordinate.
    '''
    pass

def suggest_new_residue_number_for_ligand(model, chain_id):
    '''
    Suggest a suitable residue number for a new ligand, based on what is already
    in the chain.
    '''
    from chimerax.atomic import Residue
    residues = model.residues[model.residues.chain_ids == chain_id]
    ligand_residues = residues[residues.polymer_types==Residue.PT_NONE]
    if ligand_residues:
        return max(ligand_residues.numbers)+1

    last_polymeric_residue_num = max(residues.numbers)
    # Start ligands on a nice round multiple of 1000
    ligand_num = round(last_polymeric_residue_num+1000, -3)
    if ligand_num > 9999:
        raise TypeError('The PDB format does not support residue numbers greater '
            'than 9999. Consider adding your ligand to a different chain.')
    return ligand_num
