# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 11-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll



'''
General tools for working with and manipulating OpenMM forcefields.
'''

class Template_Mgr:
    def __init__(self):
        from collections import defaultdict
        self._model_to_templates = defaultdict(dict)

    def assign_template_to_residue(self, residue, template_name):
        s = r.structure
        self.template_dict(s)[r] = template_name

    def template_dict(self, model):
        return self._model_to_template[model]

    def find_template(self, model, residue, forcefield):
        td = self.template_dict(model)
        tmpl = td.get(r, None)
        if tmpl is not None:
            # Double-check that it still matches. Set it to None if it doesn't.
            pass
        if tmpl is None:
            # Find the matching template in the forcefield definition
            pass




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

    _, sg = find_maximal_isomorphous_fragment(rn, tn)

    remainder = set(residue.atoms).difference(set(sg.keys()))
    for a in remainder:
        if a.element.name != 'H':
            raise TypeError('Residue contains at least one non-hydrogen atom not found in the template!')
    for a in remainder:
        a.delete()

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
    from networkx import isomorphism as iso
    if len(residue_graph) >= len(template_graph):
        lg = residue_graph
        sg = template_graph
        residue_larger = True
    else:
        lg = template_graph
        sg = residue_graph
        residue_larger = False

    gm = iso.GraphMatcher(lg, sg,
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
    from openmm.app import ForceField
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
    atoms = residue.atoms
    rn.add_nodes_from(atoms)
    nx.set_node_attributes(rn, {a: {'element': e} for a, e in zip(atoms, atoms.element_names)})
    rn.add_edges_from([b.atoms for b in atoms.intra_bonds])
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
    nx.set_node_attributes(tn, {a.name: {'element': a.element.symbol} for a in tatoms})
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

def heavy_atom_names_match(residue, template):
    residue_names = set(residue.atoms[residue.atoms.element_names!='H'].names)
    template_names = set(a.name for a in template.atoms if a.element.symbol != 'H')
    return (residue_names.issubset(template_names))

def is_single_fragment(residue_graph):
    import networkx as nx
    num_subgraphs = len(list(nx.connected_components(residue_graph)))
    return (num_subgraphs == 1)
