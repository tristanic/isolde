
_fetched_templates = set()

def find_incorrect_residues(session, model, heavy_atoms_only = True):
    from chimerax.atomic import Residue, mmcif, TmplResidue
    residues = model.residues
    questionable = []
    for r in residues:
        pt = r.polymer_type
        tmpl = None
        if pt in (Residue.PT_AMINO, Residue.PT_NUCLEIC):
            start=True
            end=True
            if pt == Residue.PT_AMINO:
                fa = r.find_atom('N')
                la = r.find_atom('C')
            else:
                fa = r.find_atom('P')
                la = r.find_atom("O3'")
            if fa is not None:
                for fn in fa.neighbors:
                    if fn.residue != r:
                        start=False
                        break
            if la is not None:
                for ln in la.neighbors:
                    if ln.residue != r:
                        end=False
                        break
            try:
                tmpl = TmplResidue.get_template(r.name, start=start, end=end)
            except ValueError:
                tmpl = None
        if tmpl is None:
            if r.name not in _fetched_templates:
                session.logger.info('Fetching CCD definition for residue {} {}{}'.format(
                    r.name, r.chain_id, r.number
                ))
            try:
                tmpl = mmcif.find_template_residue(session, r.name)
                _fetched_templates.add(r.name)
            except ValueError:
                session.logger.warning('Template {} not found in the Chemical Components Dictionary'.format(r.name))
                continue
        if heavy_atoms_only:
            ra_names = set(r.atoms[r.atoms.element_names != 'H'].names)
            ta_names = set([a.name for a in tmpl.atoms if a.element.name != 'H'])
        else:
            ra_names = set(r.atoms.names)
            ta_names = set([a.name for a in tmpl.atoms])
        ra_residuals = ra_names.difference(ta_names)
        ta_residuals = ta_names.difference(ra_names)
        if len(ta_residuals):
            if end:
                if pt == Residue.PT_AMINO:
                    print('C-terminal residue {} {}{}; ra_residuals: {}; ta_residuals: {}'.format(
                        r.name, r.chain_id, r.number, ra_residuals, ta_residuals
                    ))
                if pt == Residue.PT_AMINO and not len(ra_residuals) and ta_residuals == set(('OXT',)):
                    # Dangling C-terminal peptide. Allow.
                    continue
                elif pt == Residue.PT_NUCLEIC and not heavy_atoms_only and not len(ra_residuals) and ta_residuals == set(("HO5'",)):
                    # Dangling 5' end. Allow
                    continue
            if start and not heavy_atoms_only:
                if pt == Residue.PT_AMINO and ta_residuals == set(('H',)) and ra_residuals == set(('H1','H2','H3')):
                    # Dangling N-terminal peptide. Allow
                    continue
                elif pt == Residue.PT_NUCLEIC and ta_residuals == set(("HO3'",)):
                    # Dangling 3' end. Allow.
                    continue

            questionable.append(r)
    return questionable

def residue_graph(residue):
    '''
    Make a :class:`networkx.Graph` representing the connectivity of a residue's
    atoms.
    '''
    import networkx as nx
    rn = nx.Graph()
    atoms = residue.atoms
    rn.add_nodes_from(atoms)
    nx.set_node_attributes(rn, {a: {'element': e} for a, e in zip(atoms, atoms.element_names)})
    rn.add_edges_from([b.atoms for b in atoms.intra_bonds])
    return rn

def template_graph(template):
    '''
    Make a :class:`networkx.Graph` representing the connectivity of a template's
    atoms.
    '''
    import networkx as nx
    tn = nx.Graph()
    atoms = template.atoms
    tn.add_nodes_from(atoms)
    nx.set_node_attributes(tn, {a: {'element': e} for a, e in zip(atoms, [aa.element.name for aa in atoms])})
    bonds = set()
    for a in atoms:
        bonds.update(set(a.bonds))
    tn.add_edges_from([b.atoms for b in bonds])
    return tn

def find_maximal_isomorphous_fragment(residue_graph, template_graph):
    '''
    When a residue doesn't quite match its template, there can be various
    explanations:

    - the residue is incomplete (e.g. lacking hydrogens, or actually truncated)
    - the hydrogen addition has misfired, adding too many or too few hydrogens
    - this isn't actually the correct template for the residue

    This method compares the graph representations of residue and template to
    find the maximal overlap, returning a dict mapping nodes in the larger
    to their equivalent in the smaller.

    Args:
        * residue_graph
            - a :class:`networkx.Graph` object representing the residue atoms
              present in the model (e.g. as built by :func:`residue_graph`)
        * template_graph
            - a :class:`networkx.Graph` object representing the residue atoms
              present in the template (e.g. as built by :func:`template_graph`)

    Returns:
        * bool
            - True if the residue is larger than or equal to the template,
              False otherwise
        * dict
            - mapping of nodes in the larger graph (or the residue graph if the
              sizes are equal) to their equivalent in the smaller
    '''
    import networkx as nx
    from networkx import isomorphism as iso
    if len(residue_graph) >= len(template_graph):
        lg = residue_graph
        sg = template_graph
        residue_larger = True
    else:
        lg = template_graph
        sg = residue_graph
        residue_larger = False

    max_nodes = 0
    largest_sg = None
    for cg in nx.connected.connected_components(sg):
        if len(cg) > max_nodes:
            max_nodes = len(cg)
            largest_sg = sg.subgraph(cg)
    gm = iso.GraphMatcher(lg, largest_sg,
            node_match=iso.categorical_node_match('element', None))
    subgraph = None
    # GraphMatcher will typically find multiple maximal subgraphs since we're
    # filtering by element only, not by strict atom name. Since CCD and AMBER
    # don't always agree on atom nomenclature (and ChimeraX hydrogen names may
    # not necessarily match either) there is not much more we can do. Should be
    # safe to just take the first.
    for subgraph in gm.subgraph_isomorphisms_iter():
        break

    if (subgraph is None
        or len(subgraph) <2
        or len(subgraph)/len(lg) < 0.1
        ):
        raise TypeError('Residue does not appear to match template!')

    return residue_larger, subgraph

def add_metal_bonds_from_template(residue, template):
    m = residue.structure
    metal_atoms = [a for a in template.atoms if a.element.is_metal]
    for met in metal_atoms:
        rmet = residue.find_atom(met.name)
        if rmet is None:
            raise TypeError('Residue does not have the required metal ion! Use fix_residue_from_template() instead.')
        for n in met.neighbors:
            rn = residue.find_atom(n.name)
            if not rn in rmet.neighbors:
                m.new_bond(rmet, rn)

def fix_residue_from_template(residue, template):
    m = residue.structure
    from chimerax.atomic.struct_edit import add_dihedral_atom
    from chimerax.core.geometry import distance, angle, dihedral
    residue.atoms[residue.atoms.element_names=='H'].delete()
    rnames = set(residue.atoms.names)
    tnames = set([a.name for a in template.atoms])
    wrong_names = rnames.difference(tnames)
    for wn in wrong_names:
        residue.find_atom(wn).delete()
        rnames.remove(wn)

    # Add missing bonds between atoms already in the residue
    for ta in template.atoms:
        ra = residue.find_atom(ta.name)
        if ra is not None:
            for tn in ta.neighbors:
                rn = residue.find_atom(tn.name)
                if rn is not None and rn not in ra.neighbors:
                    m.new_bond(ra, rn)

    rg = residue_graph(residue)
    tg = template_graph(template)
    residue_larger, matched_nodes = find_maximal_isomorphous_fragment(rg, tg)

    if len(matched_nodes) < 3:
        raise TypeError('Need at least 3 contiguous atoms to complete residue!')

    # Delete any isolated atoms and rebuild from template
    from chimerax.atomic import Atoms
    if residue_larger:
        conn_ratoms = Atoms(matched_nodes.keys())
    else:
        conn_ratoms = Atoms(matched_nodes.values())
    residue.atoms.subtract(conn_ratoms).delete()
    built = set(conn_ratoms)
    all_tatom_names = set([a.name for a in template.atoms])
    while len(residue.atoms) < len(template.atoms):
        tnames = [a.name for a in matched_nodes.keys()]
        rnames = [a.name for a in matched_nodes.values()]
        print('Atom numbers: Residue {}, Template {}'.format(
            len(rnames), len(tnames)
        ))
        m.session.logger.status('Still missing: {}'.format(
            ','.join(all_tatom_names.difference(set(rnames)))
        ))
        new_atoms = {}
        for ta, ra in matched_nodes.items():
            for tn in ta.neighbors:
                if not residue.find_atom(tn.name):
                    a = build_next_atom_from(tn.name, ta.name, residue, template)
                    new_atoms[tn] = a
        matched_nodes = new_atoms
    bonds = set()
    for a in template.atoms:
        bonds.update(set(a.bonds))
    for b in bonds:
        a1, a2 = [residue.find_atom(a.name) for a in b.atoms]
        if a2 not in a1.neighbors:
            m.new_bond(a1, a2)


def build_next_atom_from(next_atom_name, stub_atom_name, residue, template):
    from chimerax.atomic import struct_edit
    from chimerax.core.geometry import distance, angle, dihedral
    r = residue
    m = r.structure
    tnext = template.find_atom(next_atom_name)
    if tnext is None:
        raise TypeError('Template does not contain an atom with that name!')
    tstub = template.find_atom(stub_atom_name)
    rstub = r.find_atom(stub_atom_name)
    # paired_atom_count = 1
    # paired_atoms = {tstub: rstub,}
    # search_list = tstub.neighbors
    # while paired_atom_count < 4:
    #     seen = set()
    #     new_search_list = set()
    #     for ta in search_list:
    #         if ta.element.name != 'H' and ta not in seen:
    #             seen.add(ta)
    #             new_search_list.update(set(ta.neighbors))
    #             ra = r.find_atom(ta.name)
    #             if ra is not None:
    #                 paired_atoms[ta] = ra
    #                 paired_atom_count += 1
    #     search_list = new_search_list

    n1 = rstub
    n2 = n3 = None
    t_direct_neighbors = []
    r_direct_neighbors = []
    for a2 in tstub.neighbors:
        if a2.element.name != 'H':
            n2 = r.find_atom(a2.name)
            if n2:
                t_direct_neighbors.append(a2)
                r_direct_neighbors.append(n2)
    if len(t_direct_neighbors) > 1:
        a2, a3 = t_direct_neighbors[:2]
        n2, n3 = r_direct_neighbors[:2]
    else:
        a2 = t_direct_neighbors[0]
        n2 = r_direct_neighbors[0]
    if not n2:
        raise TypeError('No n2 found - Not enough connected atoms to form a dihedral!')
    if not n3:
        for a3 in a2.neighbors:
            if a3 not in (a2, tstub) and a3.element.name != 'H':
                n3 = r.find_atom(a3.name)
                if n3:
                    break
        if not n3:
            raise TypeError('No n3 found - Not enough connected atoms to form a dihedral!')
    # from chimerax.core.geometry import align_points
    # import numpy

    #tas, ras = paired_atoms.
    # tf = align_points(numpy.array([a.coord for a in paired_atoms.keys()]), numpy.array([a.coord for a in paired_atoms.values()]))[0]

    # tf = align_points(numpy.array([a.coord for a in [tstub, a2, a3]]),
    #     numpy.array([a.coord for a in [rstub, n2, n3]]))[0]
    # apos = tf*tnext.coord
    # a = struct_edit.add_atom(next_atom_name, tnext.element, residue, apos)
    # return a


    dist = distance(tnext.coord, tstub.coord)
    ang = angle(tnext.coord, tstub.coord, a2.coord)
    dihe = dihedral(tnext.coord, tstub.coord, a2.coord, a3.coord)
    # print('{}: {} {} {}'.format(next_atom_name, dist, ang, dihe))
    a = struct_edit.add_dihedral_atom(next_atom_name, tnext.element, n1, n2, n3, dist, ang, dihe)
    return a
