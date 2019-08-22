
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
