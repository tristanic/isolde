def simple_symmetry(session, model, reference_chain_ids, symmetry_list):
    import numpy
    from chimerax.core.commands import align
    from chimerax.core.atomic import Residues, Atoms, concatenate
    def atoms_sorted_by_chain(chains):
        sort_indices = numpy.argsort(chains.chain_ids)
        ret = Atoms()
        for i in sort_indices:
            r = Residues(chains[i].residues)
            ret = concatenate([ret, r.atoms])
        return ret
    m = model
    refchains = m.chains[numpy.in1d(m.chains.chain_ids, reference_chain_ids)]
    ref = atoms_sorted_by_chain(refchains)
    ref_pos = ref.coords
    for chains in symmetry_list:
        symchains = m.chains[numpy.in1d(m.chains.chain_ids, chains)]
        sym = atoms_sorted_by_chain(symchains)
        align.align(session, ref, to_atoms = sym, move = ref)
        sym.coords = ref.coords
    ref.coords = ref_pos

def simple_symmetry(session, model, reference_chain_ids, symmetry_list):
    import numpy
    from chimerax.core.commands import align
    from chimerax.core.atomic import Residues, Atoms, concatenate
    def atoms_sorted_by_chain(chains):
        sort_indices = numpy.argsort(chains.chain_ids)
        ret = Atoms()
        for i in sort_indices:
            r = Residues(chains[i].residues)
            ret = concatenate([ret, r.atoms])
        return ret
    m = model
    refchains = m.chains[numpy.in1d(m.chains.chain_ids, reference_chain_ids)]
    ref = atoms_sorted_by_chain(refchains)
    ref_pos = ref.coords
    for chains in symmetry_list:
        symchains = m.chains[numpy.in1d(m.chains.chain_ids, chains)]
        sym = atoms_sorted_by_chain(symchains)
        align.align(session, ref, to_atoms = sym, move = ref)
        sym.coords = ref.coords
    ref.coords = ref_pos
