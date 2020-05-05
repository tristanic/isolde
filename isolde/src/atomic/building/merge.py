from chimerax.isolde.atomic.building import set_new_atom_style

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
    if (chain_id or anchor_n or anchor_c or (renumber_from is not None)) \
            and len(residues.unique_chain_ids) != 1:
        raise UserError('If reassigning chain ID, renumbering or specifying '
            'N- and/or C-terminal anchors, all residues to be copied must be '
            'from a single chain!')
    from chimerax.atomic import Residue, Residues, Atoms
    if (anchor_n or anchor_c):
        protein_residues = residues[residues.polymer_types==Residue.PT_AMINO]
        if not len(protein_residues):
            raise UserError('N- and/or C-terminal anchors were specified, but '
                'the copied selection does not contain any amino acid residues!')

    import numpy
    fm = us[0]
    m = target_model
    tpbg = m.pseudobond_group('missing structure')
    residues = Residues(sorted(residues, key=lambda r: (r.chain_id, r.number, r.insertion_code)))
    atoms = residues.atoms
    coords = atoms.coords
    atom_map = {}
    if renumber_from:
        offset = residues[0].number - renumber_from
    else:
        offset = 0

    def new_residue_number(r):
        if r in residues:
            return r.number+offset
        return r.number
    merged_residues = m.residues.merge(residues)
    merged_residues = Residues(sorted(merged_residues,
        key=lambda r: (r.chain_id, new_residue_number(r), r.insertion_code)
        ))
    new_residue_mask = numpy.in1d(merged_residues, residues)
    new_residue_indices = numpy.argwhere(new_residue_mask).ravel()
    existing_residue_mask = numpy.in1d(merged_residues, m.residues)
    existing_residue_indices = numpy.argwhere(existing_residue_mask).ravel()

    insertion_point_map = {}

    if chain_id is not None:
        cids = [chain_id]
    else:
        cids = residues.unique_chain_ids
    for cid in cids:
        existing_residues = m.residues[m.residues.chain_ids==cid]
        if not len(existing_residues):
            insertion_point_map[cid] = None
            continue
        existing_residue_numbers = numpy.array([str(r.number)+r.insertion_code for r in existing_residues])
        cres = residues[residues.chain_ids==cid]
        new_residue_numbers = numpy.array([str(r.number+offset)+r.insertion_code for r in cres])

        duplicate_flags = numpy.in1d(new_residue_numbers, existing_residue_numbers)
        if numpy.any(duplicate_flags):
            dup_residues = cres[duplicate_flags]
            err_str = ('The requested merge could not be completed because the '
                'following residues in chain {} (after applying any renumbering) '
                'will have the same residue numbers as existing residues in '
                'the target: {}'
            ).format(cid, ', '.join(str(r.number)+r.insertion_code for r in dup_residues))
            raise UserError(err_str)

        chain_mask = merged_residues.chain_ids == cid
        new_r_in_chain_mask = numpy.logical_and(chain_mask, new_residue_mask)
        last_new_r_index = numpy.argwhere(new_r_in_chain_mask)[-1]
        greater_indices = existing_residue_indices[existing_residue_indices > last_new_r_index]
        if len(greater_indices):
            insertion_point_map[cid] = merged_residues[greater_indices[0]]
        else:
            insertion_point_map[cid] = None


    current_residue = None

    first_index = new_residue_indices[0]
    if first_index > 0:
        prev_res = merged_residues[first_index-1]
    else:
        prev_res = None
    prev_new_res = None

    for merged_index, r in zip(new_residue_indices, residues):
        if chain_id:
            cid = chain_id
        else:
            cid = r.chain_id
        precedes = insertion_point_map[cid]
        insertion_code = r.insertion_code
        if insertion_code=='':
            insertion_code = ' '
        nr = m.new_residue(r.name, cid, r.number-offset, insert=insertion_code,
            precedes=precedes)
        nr.ribbon_hide_backbone = r.ribbon_hide_backbone
        nr.ribbon_display = r.ribbon_display
        nr.ribbon_color = r.ribbon_color
        nr.ss_type = r.ss_type

        for a in r.atoms:
            na = atom_map[a] = m.new_atom(a.name, a.element)
            na.coord = a.coord
            na.bfactor = a.bfactor
            na.aniso_u6 = a.aniso_u6
            na.occupancy = a.occupancy
            na.draw_mode = a.draw_mode
            na.color = a.color
            na.display = a.display
            nr.add_atom(na)
            for n in a.neighbors:
                nn = atom_map.get(n, None)
                if nn is not None:
                    m.new_bond(na, nn)


        if prev_res is not None:
            if (
              r.polymer_type==Residue.PT_AMINO and prev_res not in r.neighbors
              and prev_res.chain_id == cid
              and prev_res.polymer_type == Residue.PT_AMINO):
                if prev_res.structure != r.structure:
                    if precedes.chain_id == cid and precedes.polymer_type == Residue.PT_AMINO:
                        ratoms = prev_res.atoms.merge(precedes.atoms)
                        tpbg.pseudobonds[tpbg.pseudobonds.between_atoms(ratoms)].delete()
                    pc = prev_res.find_atom('C')
                else:
                    if prev_new_res is not None:
                        pc = prev_new_res.find_atom('C')
                    else:
                        pc = None
                nn = nr.find_atom('N')
                if pc and nn:
                    tpbg.new_pseudobond(pc, nn)
                if precedes.polymer_type==Residue.PT_AMINO:
                    nc = nr.find_atom('C')
                    pn = precedes.find_atom('N')
                    if nc and pn:
                        tpbg.new_pseudobond(nc, pn)
        prev_res = r
        prev_new_res = nr
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
    if update_style:
        new_atoms.displays=True
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
