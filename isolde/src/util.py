# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 27-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



import numpy
from chimerax.atomic import Residue, AtomicStructure
def is_continuous_protein_chain(sel, allow_single = False):
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
    if len(indices) <= 1 and not allow_single:
        return False

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
        p = p[0]
        first_index = p.index(r0)
        if first_index != -1:
            return p
    raise IndexError('Polymer not found!')


def add_disulfides_from_model_metadata(model):
    m_id = model.id_string
    from chimerax.core.commands import atomspec
    metadata = model.metadata
    try:
        disulfide_list = metadata['SSBOND']
    except KeyError:
        return
    for disulfide in disulfide_list:
        sym1 = None
        sym2 = None
        d = disulfide.split()
        chain1, res1 = d[3], d[4]
        chain2, res2 = d[6], d[7]
        if len(d) > 8:
            sym1, sym2 = d[8], d[9]
        if sym1 is not None and sym1 != sym2:
            # Disulfide across a symmetry interface. Ignore for now.
            continue
        arg = atomspec.AtomSpecArg('ss')
        thearg = '#{}/{}:{}@{}|#{}/{}:{}@{}'.format(
            m_id, chain1, res1, 'SG', m_id, chain2, res2, 'SG')
        aspec = arg.parse(thearg, model.session)
        atoms = aspec[0].evaluate(model.session).atoms
        bonds = atoms.intra_bonds
        if not len(bonds):
            a1, a2 = atoms
            b = model.new_bond(a1, a2)

def compiled_lib_extension():
    import platform
    pname = platform.system()
    if pname == "Windows":
        return "dll"
    elif pname == "Darwin":
        return "dylib"
    return "so"
