from . import set_new_atom_style

def merge_fragment(target_model, residues, chain_id=None, renumber_from=None,
        anchor_n=None, anchor_c=None, update_style=True, transform=None):
    '''
    Copy the atoms from a fragment into the current model, optionally reassigning
    chain ID and numbers. All residues in the fragment must have the same chain
    ID. If alternate conformers are present, only the active one will be copied.

    * Args:

        - target_model: the model to copy into
        - residues: the residues to be copied. Isotropic B-factors will be
            copied, but aniso_u records will be ignored.
        - chain_id: if provided, all residues must be from one chain. The copied
            residues will be given this chain ID.
        - renumber_from: if provided, all residues will be renumbered with an
            offset of (lowest residue number - renumber_from)
        - anchor_n: an amino acid residue or None. If provided, the first amino
            acid residue in the fragment will be linked to this one. Throws an
            error if anchor_n has another residue linked at its C-terminus.
        - anchor_c: an amino acid residue or None. If provided, the last amino
            acid residue in the fragment will be linked to this one. Throws an
            error if anchor_c has another residue linked at its C-terminus.
        - transform: a Place or None. If provided, the atoms will be placed at
            the transformed coordinates.
    '''
    from chimerax.core.errors import UserError
    us = residues.unique_structures
    if len(us) != 1:
        raise UserError('All residues to be copied must be from the same model!')
    fm = us[0]
    fpbg = fm.pseudobond_group('missing structure')
    m = target_model
    tpbg = m.pseudobond_group('missing structure')
    if (chain_id or anchor_n or anchor_c or (renumber_from is not None)) \
            and len(residues.unique_chain_ids) != 1:
        raise UserError('If reassigning chain ID, renumbering or specifying '
            'N- and/or C-terminal anchors, all residues to be copied must be '
            'from a single chain!')
    from chimerax.atomic import Residue, Residues, Atoms
    residues = Residues(sorted(residues, key=lambda r: r.number))
    if (anchor_n or anchor_c):
        from chimerax.atomic import Residue
        protein_residues = residues[residues.polymer_types==Residue.PT_AMINO]
        if not len(protein_residues):
            raise UserError('N- and/or C-terminal anchors were specified, but '
                'the copied selection does not contain any amino acid residues!')
    atoms = residues.atoms
    coords = atoms.coords
    atom_map = {}
    if renumber_from:
        offset = residues[0].number - renumber_from
    else:
        offset = 0
    current_residue = None
    cid = chain_id
    for a, coord in zip(atoms, coords):
        if a.residue != current_residue:
            r = a.residue
            current_residue = r
            if chain_id:
                cid = chain_id
            else:
                cid = r.chain_id
            insertion_code = r.insertion_code
            if insertion_code=='':
                insertion_code = ' '
            nr = m.new_residue(r.name, cid, r.number-offset, insert=insertion_code)
            nr.ss_type = r.ss_type
        na = atom_map[a] = m.new_atom(a.name, a.element)
        na.coord = coord
        na.bfactor = a.bfactor
        na.aniso_u6 = a.aniso_u6
        na.occupancy = a.occupancy
        nr.add_atom(na)
        for n in a.neighbors:
            nn = atom_map.get(n, None)
            if nn is not None:
                m.new_bond(na, nn)
    new_atoms = Atoms(list(atom_map.values()))
    if transform is not None:
        # Using Atoms.transform() rather than simply transforming the coords,
        # because this also correctly transforms any anisotropic B-factors.
        new_atoms.transform(transform)

    if anchor_n:
        anchor_atom = anchor_n.find_atom('C')
        link_atom = atom_map[protein_residues[0].find_atom('N')]
        _remove_excess_terminal_atoms(anchor_atom)
        _remove_excess_terminal_atoms(link_atom)
        m.new_bond(anchor_atom, link_atom)
    if anchor_c:
        anchor_atom = anchor_c.find_atom('N')
        link_atom = atom_map[protein_residues[-1].find_atom('C')]
        _remove_excess_terminal_atoms(anchor_atom)
        _remove_excess_terminal_atoms(link_atom)
        m.new_bond(anchor_atom, link_atom)
    f_pbonds = fpbg.pseudobonds
    for pb in f_pbonds:
        pb_atoms = [atom_map.get(a) for a in pb.atoms]
        if None not in pb_atoms:
            tpbg.new_pseudobond(*pb_atoms)
    if update_style:
        set_new_atom_style(m.session, new_atoms)
    return new_atoms

def _remove_excess_terminal_atoms(atom):
    if atom.name=='C':
        for n in atom.neighbors:
            if n.name=='OXT':
                n.delete()
    elif atom.name=='N':
        for n in atom.neighbors:
            if n.name in ('H2', 'H3'):
                n.delete()
            elif n.name == 'H1':
                n.name = 'H'
