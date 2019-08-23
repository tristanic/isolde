
def create_ncs_copy(model, residues, transform, new_chain_id = None):
    '''
    Create a NCS copy of a chain using the given rigid-body transform. Residues
    must all be from a single chain. If no chain ID is specified, the next
    available chain ID will be used.
    '''
    if new_chain_id is None:
        from .build_utils import next_chain_id
        new_chain_id = next_chain_id(model)

    from .merge import merge_fragment
    return merge_fragment(model, residues, chain_id=new_chain_id, transform=transform)

def create_all_ncs_copies(model, transforms):
    '''
    Expand the model to include all NCS copies.
    '''
    chain_ids = model.chains.chain_ids
    residues = model.residues
    chains = [residues[residues.chain_ids==cid] for cid in chain_ids]
    for tf in transforms:
        if not tf.is_identity():
            for chain in chains:
                create_ncs_copy(model, chain, tf)

def apply_strict_ncs(template_atoms, ncs_map):
    '''
    Copy the transformed coordinates of the template atoms to the given NCS
    atoms.

    Args:

        * template_atoms: a :class:`Atoms` instance
        * ncs_map: a list of (:class:`Place`, :class:`Atoms`) pairs. Each
            array of atoms must have a 1:1 correspondence to the template_atoms.
            Atoms will be sorted according to (chain, residue number, insertion
            code, atom name) before applying the NCS operators.
    '''
    from chimerax.atomic import Atoms
    template_atoms = Atoms(sorted(template_atoms, key=lambda a:
        (a.residue.chain_id, a.residue.number, a.residue.insertion_code, a.name)))
    for (place, ncs_atoms) in ncs_map:
        ncs_atoms = Atoms(sorted(ncs_atoms, key=lambda a:
            (a.residue.chain_id, a.residue.number, a.residue.insertion_code, a.name)))
        ncs_atoms.coords = place*template_atoms.coords
