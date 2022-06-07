# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 04-Jan-2021
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



import numpy
from chimerax.atomic import Residue, AtomicStructure
from contextlib import contextmanager

@contextmanager
def block_trigger_handler(handler):
    name = handler._name
    func = handler._func
    triggerset = handler._trigger_set
    triggerset.remove_handler(handler)
    try:
        yield
    finally:
        triggerset.add_handler(name, func)

@contextmanager
def block_managed_trigger_handler(object, handler_name):
    '''
    `object.handler_name` should be a ChimeraX trigger handler. The original handler 
    will be deleted, and replaced with an identical new one on exit of the context.
    '''
    h = getattr(object, handler_name)
    name = h._name
    func = h._func
    triggerset = h._trigger_set
    triggerset.remove_handler(h)
    try:
        yield
    finally:
        setattr(object, handler_name, triggerset.add_handler(name, func))


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

def peptide_has_adducts(residues):
    from chimerax.atomic import Residue
    residues = residues[residues.polymer_types==Residue.PT_AMINO]
    for r in residues:
        neighbors = r.neighbors
        if len(neighbors) > 2:
            return True
        for n in neighbors:
            bonds = n.bonds_between(r)
            for b in bonds:
                for a in b.atoms:
                    if a.residue==r and a.name not in ('C', 'N'):
                        return True
    return False



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
    from chimerax.atomic.struct_edit import add_bond
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
            add_bond(a1, a2)

def expand_selection(residues, num_steps):
    from .molobject import c_function
    import ctypes
    f = c_function('expand_selection',
        args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int32)
        )
    f(residues._c_pointers, len(residues), num_steps)

def contract_selection(residues, num_steps):
    from chimerax.atomic import Residues
    for _ in range(num_steps):
        removed = []
        for r in residues:
            selected_neighbors = 0
            for n in r.neighbors:
                if n.selected:
                    selected_neighbors += 1
            if selected_neighbors == 1:
                removed.append(r)
        removed = Residues(removed)
        removed.atoms.selected=False
        removed.atoms.intra_bonds.selected=False
        for r in removed:
            for n in r.neighbors:
                r.bonds_between(n).selected=False
        residues = residues.subtract(removed)
                


def compiled_lib_extension():
    import platform
    pname = platform.system()
    if pname == "Windows":
        return "dll"
    elif pname == "Darwin":
        return "dylib"
    return "so"

class SafeTempDir:
    '''
    Create and change to a temporary directory. When called with

    with SafeTempDir():
        do_something()

    ... it is guaranteed to change back to the original working directory and
    delete the temporary directory.
    '''
    def __init__(self):
        import os
        self.cwd = os.path.abspath(os.curdir)

    def __enter__(self):
        import tempfile, os
        td = self._temp_dir = tempfile.mkdtemp()
        os.chdir(td)

    def __exit__(self, exc_type, exc_value, exc_traceback):
        import os
        os.chdir(self.cwd)
        import shutil
        shutil.rmtree(self._temp_dir)
