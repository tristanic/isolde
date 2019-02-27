'''
General tools for working with and manipulating OpenMM forcefields.
'''

class Template_Mgr:
    def __init__(self):
        from collections import defaultdict
        self._model_to_templates = defaultdict(dict)

    def assign_template_to_residue(self, residue, template_name):
        s = r.structure
        self._model_to_template[s][r] = template_name

    def template_dict(self, model):
        return self._model_to_template[model]


def delete_extraneous_hydrogens(residue, template):
    '''
    Addition of hydrogens to a model is still something of a black art,
    particularly when attempting to guess the chemistry directly from the
    geometry of the heavy atoms. Outside of the standard protein/nucleic acid
    framework, tools such as ChimeraX's `AddH` occasionally get things wrong.
    The most common issue is too many hydrogens, which is relatively simple to
    fix.
    '''
    if len(residue.atoms) == len(template.atoms):
        return
    rn = residue_graph(residue)
    tn = template_graph(template)

    import networkx as nx
    gm = nx.isomorphism.GraphMatcher(rn, tn)
    sg = None
    for sg in gm.subgraph_isomorphisms_iter():
        break

    if not sg:
        raise TypeError('Template {} does not appear to match residue {}'.format(
            template.name, residue.name
        ))

    remainder = set(residue.atoms).difference(set(sg.keys()))
    for a in remainder:
        if a.element.name != 'H':
            raise TypeError('Residue contains at least one non-hydrogen atom not found in the template!')
    for a in remainder:
        a.delete()


def make_truncated_residue_template(residue, full_template):
    '''
    Given a :class:`chimerax.Residue` representing an incompletely-built residue
    :class:`simtk.openmm.app.forcefield.ForceField._TemplateData` template,
    create a new partial template covering only the atoms present. The resulting
    template will of course be highly artificial, and should *only* be used in
    the context of model-building in medium-high resolution maps. At very low
    resolution, it is almost always better to keep all residues complete.

    Due to the artificial nature of this approach, some strict rules apply:

        - the partial residue must be a single, contiguous fragment
        - any rings must be either complete, or completely missing
        - the names of all heavy atoms in the residue must match their
          counterparts in the template
    '''
    rn = residue_graph(residue)
    tn = template_graph(full_template)

    if not is_single_fragment(rn):
        raise TypeError('Partial residue must be a single fragment!')
    if not heavy_atom_names_match(residue, full_template):
        raise TypeError('To model a partial residue, the names of all heavy '
            'atoms must match their counterparts in the template!')
    if any_incomplete_rings(residue, tn):
        raise TypeError('Cannot model a partial ring!')

    keep_heavy_atom_names = set(residue.atoms[residue.atoms.element_names != 'H'].names)
    from simtk.openmm.app import ForceField
    new_template = ForceField._TemplateData('TRUNC_'+residue.name+residue.chain_id+str(residue.number))
    for atom in full_template.atoms:
        if atom.name in keep_heavy_atom_names:
            new_template.addAtom(ForceField._TemplateAtomData(atom.name, atom.type, atom.element))
            new_template.atoms[-1].externalBonds = atom.externalBonds
    # Find the hydrogens to copy over
    for (a1, a2) in full_template.bonds:
        atom1, atom2 = full_template.atoms[a1], full_template.atoms[a2]
        names = set((atom1.name, atom2.name))
        if not len(names.intersection(keep_heavy_atom_names))==1:
            continue
        for atom in (atom1, atom2):
            if atom.element.symbol == 'H':
                new_template.addAtom(ForceField._TemplateAtomData(atom.name, atom.type, atom.element))
    new_template_atom_names = set(a.name for a in new_template.atoms)
    for (a1, a2) in full_template.bonds:
        names = set((full_template.atoms[a1].name, full_template.atoms[a2].name))
        if names.issubset(new_template_atom_names):
            new_template.addBondByName(*names)

    # AddH will almost invariably add extra hydrogens to the truncation point in
    # truncated residues, so may as well delete those here
    delete_extraneous_hydrogens(residue, new_template)
    return new_template

def residue_graph(residue):
    '''
    Make a :class:`networkx.Graph` representing the connectivity of a residue's
    atoms.
    '''
    import networkx as nx
    rn = nx.Graph()
    rn.add_nodes_from(residue.atoms)
    rn.add_edges_from([b.atoms for b in residue.atoms.intra_bonds])
    return rn

def template_graph(template):
    '''
    Make a :class:`networkx.Graph` representing the connectivity of an OpenMM
    template.
    '''
    import networkx as nx
    tn = nx.Graph()
    tatoms = template.atoms
    tn.add_nodes_from([a.name for a in tatoms])
    tn.add_edges_from([(tatoms[a1].name, tatoms[a2].name) for (a1, a2) in template.bonds])
    return tn

def is_valid_partial_fragment(residue, template):
    '''
    Check if it is possible to simulate a partially-modelled residue by
    truncating its forcefield template
    '''
    return (is_single_fragment(residue) and heavy_atoms_match(residue, template)
        and not any_incomplete_rings(residue, template))

def any_incomplete_rings(residue, template_graph):
    atom_names = set(residue.atoms.names)
    rings = find_cycles(template_graph)
    for ring in rings:
        ratoms = atom_names.intersection(ring)
        if len(ratoms) and not ratoms == ring:
            return True
    return False


def find_cycles(template_graph):
    import networkx as nx
    return [set(alist) for alist in nx.cycle_basis(template_graph)]

# def find_cycles(residue):
#     import networkx as nx
#     g = nx.Graph()
#     g.add_nodes_from(residue.atoms)
#     g.add_edges_from([b.atoms for b in residue.atoms.intra_bonds])
#     return [set(alist) for alist in nx.cycle_basis(g)]

def heavy_atom_names_match(residue, template):
    residue_names = set(residue.atoms[residue.atoms.element_names!='H'].names)
    template_names = set(a.name for a in template.atoms if a.element.symbol != 'H')
    return (residue_names.issubset(template_names))


def is_single_fragment(residue_graph):
    import nx
    num_subgraphs = len(list(nx.connected_components(residue_graph)))
    return (num_subgraphs == 1)

    #
    # all_atoms = set(residue.atoms)
    # traversed_atoms = set()
    # _traverse_residue(residue.atoms[0], all_atoms, traversed_atoms)
    # # return (traversed_atoms == all_atoms)
    # if traversed_atoms != all_atoms:
    #     return False
    # return True

# def _traverse_residue(atom, all_atoms, traversed_atoms):
#     traversed_atoms.add(atom)
#     for a in atom.neighbors:
#         if a not in all_atoms or a in traversed_atoms:
#             continue
#         _traverse_residue(a, all_atoms, traversed_atoms)
