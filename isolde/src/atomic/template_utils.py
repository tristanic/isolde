
from chimerax.core.errors import UserError

_fetched_templates = set()

def find_incorrect_residues(session, model, heavy_atoms_only = True):
    from chimerax.atomic import Residue, TmplResidue
    from chimerax import mmcif
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

class AtomProxy:
    def __init__(self, atom):
        self.atom = atom
    def __lt__(self, other):
        return self.atom.name < other.atom.name

def _nx_residue_graph(residue):
    '''
    Make a :class:`networkx.Graph` representing the connectivity of a residue's
    atoms.
    '''
    import networkx as nx
    rn = nx.Graph()
    atoms = residue.atoms
    proxies = [AtomProxy(atom) for atom in atoms]
    proxy_map = {atom: proxy for atom, proxy in zip(atoms, proxies)}
    rn.add_nodes_from(proxies)
    nx.set_node_attributes(rn, {a: {'element': e, 'name': n} for a, e, n in zip(proxies, atoms.element_names, atoms.names)})
    rn.add_edges_from([[proxy_map[a] for a in b.atoms] for b in atoms.intra_bonds])
    return rn

def residue_graph(residue, label='element'):
    from chimerax.isolde.graph import make_graph_from_residue
    return make_graph_from_residue(residue, label=label)

def _nx_template_graph(template):
    '''
    Make a :class:`networkx.Graph` representing the connectivity of a template's
    atoms.
    '''
    import networkx as nx
    tn = nx.Graph()
    atoms = template.atoms
    proxies = [AtomProxy(atom) for atom in atoms]
    proxy_map = {atom: proxy for atom, proxy in zip(atoms, proxies)}
    tn.add_nodes_from(proxies)
    nx.set_node_attributes(tn, {a: {'element': e, 'name': n} for a, e, n in zip(proxies, [aa.element.name for aa in atoms], [aa.name for aa in atoms])})
    bonds = set()
    for a in atoms:
        bonds.update(set(a.bonds))
    tn.add_edges_from([[proxy_map[a] for a in b.atoms] for b in bonds])
    return tn

def template_graph(template, label='element'):
    from chimerax.isolde.graph import make_graph_from_residue_template
    return make_graph_from_residue_template(template, label=label)

def find_maximal_isomorphous_fragment(residue, template, match_by='element',
        limit_template_indices=None):
    '''
    When a residue doesn't quite match its template, there can be various
    explanations:

    - the residue is incomplete (e.g. lacking hydrogens, or actually truncated)
    - the hydrogen addition has misfired, adding too many or too few hydrogens
    - this isn't actually the correct template for the residue

    This method compares the graph representations of residue and template to
    find the maximal overlap, returning a dict mapping atoms in the residue to
    those in the template.

    Args:
        * residue
            - a :class:`chimerax.atomic.Residue` object
        * template
            - a :class:`chimerax.atomic.cytmpl.TmplResidue` object

    Returns:
        * matched_atoms (dict)
            - mapping of residue atoms to matching template atoms
        * residue_extra_atoms (list)
            - atoms that are in the residue but were not matched to the template
              (if the residue is disconnected into two or more fragments, all
              atoms that aren't in the largest will be found here).
        * template_extra_atoms (list)
            - atoms in the template that aren't in the residue.
    '''
    rg = residue_graph(residue, label=match_by)
    tg = template_graph(template, label=match_by)
    residue_indices, template_indices, timed_out = rg.maximum_common_subgraph(tg, big_first=True)
    if timed_out:
        ri2, ti2, to2 = rg.maximum_common_subgraph(tg, timeout=5)
        if len(ri2) > len(residue_indices):
            residue_indices, template_indices = ri2, ti2
        if to2:
            residue.session.logger.warning('Timed out trying to match residue {} to template {}. Match may not be ideal.'.format(residue.name, template.name))
    tatoms = template.atoms
    if limit_template_indices is not None:
        import numpy
        tatoms = [tatoms[i] for i in range(len(tatoms)) if i in limit_template_indices]
        mask = numpy.in1d(template_indices, limit_template_indices)
        residue_indices = residue_indices[mask]
        template_indices = template_indices[mask]

    amap = {residue.atoms[ri]: template.atoms[ti] for ri, ti in zip(residue_indices, template_indices)}

    from chimerax.atomic import Atoms
    residue_extra_atoms = Atoms(set(residue.atoms).difference(set(amap.keys())))
    template_extra_atoms = list(set(tatoms).difference(set(amap.values())))



    return amap, residue_extra_atoms, template_extra_atoms

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

def fix_residue_from_template(residue, template, rename_atoms_only=False,
        rename_residue=False, match_by='name', template_indices=None):
    import numpy
    from chimerax.atomic import Atoms
    if any([numpy.any(numpy.isnan(a.coord)) for a in template.atoms]):
        raise TypeError('Template is missing one or more atom coordinates!')
    matched_nodes, residue_extra, template_extra = find_maximal_isomorphous_fragment(residue, template, limit_template_indices=template_indices, match_by=match_by)

    if len(matched_nodes) < 3:
        raise UserError('Residue {} {}{} has only {} connected atoms in common with template {}. At least 3 matching atoms are needed.'.format(
            residue.name, residue.chain_id, residue.number, len(matched_nodes), template.name
        ))

    m = residue.structure
    session = m.session
    from chimerax.atomic.struct_edit import add_dihedral_atom
    from chimerax.geometry import distance, angle, dihedral
    # if not rename_atoms_only:
    #     residue.atoms[residue.atoms.element_names=='H'].delete()
    rnames = set(residue.atoms.names)
    tnames = set([a.name for a in template.atoms])
    # wrong_names = rnames.difference(tnames)
    # for wn in wrong_names:
    #     residue.find_atom(wn).delete()
    #     rnames.remove(wn)

    # # Add missing bonds between atoms already in the residue
    # if not rename_atoms_only:
    #     for ta in template.atoms:
    #         ra = residue.find_atom(ta.name)
    #         if ra is not None:
    #             for tn in ta.neighbors:
    #                 rn = residue.find_atom(tn.name)
    #                 if rn is not None and rn not in ra.neighbors:
    #                     m.new_bond(ra, rn)



    # Delete any isolated atoms and rebuild from template
    if len(residue_extra):
        if rename_atoms_only:
            session.logger.warning('The following atoms in {} {}{} did not match template {}, and will not be renamed: {}.'.format(
                residue.name, residue.chain_id, residue.number, template.name,
                ', '.join(residue_extra.names)
            ))
        else:
            session.logger.info('Deleted the following atoms from residue {} {}{}{}: {}'.format(
                residue.name, residue.chain_id, residue.number, residue.insertion_code, ', '.join(residue_extra.names)
            ))
            residue_extra.delete()

    conn_ratoms = Atoms(matched_nodes.keys())
    renamed_atoms = []
    for ratom, tatom in matched_nodes.items():
        if ratom.name != tatom.name:
            renamed_atoms.append((ratom.name, tatom.name))
            ratom.name = tatom.name
    if len(renamed_atoms):
        warn_str = '{} atoms were automatically renamed to match the template: '.format(len(renamed_atoms))
        warn_str += ', '.join(['->'.join(apair) for apair in renamed_atoms])
        session.logger.warning(warn_str)
    if rename_atoms_only:
        return

    built = set(conn_ratoms)
    all_tatom_names = set([a.name for a in template.atoms])
    still_missing = set(template_extra)
    remaining = len(still_missing)
    # last_len = len(residue.atoms)-1
    found = set()
    while remaining:
        for ta in still_missing:
            found_neighbors = []
            for tn in ta.neighbors:
                ra = residue.find_atom(tn.name)
                if ra:
                    found_neighbors.append((ra, tn))
            if len(found_neighbors) == 0:
                continue
            else:
                ra, tn = found_neighbors[0]
                build_next_atom_from_geometry(residue, ra, tn, ta)
            # else:
            #     build_next_atom_from_coords(residue, found_neighbors, ta)
            found.add(ta)
        still_missing = still_missing.difference(found)
        if len(still_missing) and not len(found):
            raise RuntimeError('Failed to add atoms on last iteration for residue {}. Still missing: {}'.format(
                '{}{}'.format(residue.chain_id, residue.number), ','.join(still_missing)
            ))
        # m.session.logger.status('Still missing: {}'.format(
        #     ','.join([ta.name for ta in still_missing])
        # ))
        remaining = len(still_missing)
    bonds = set()
    for a in template.atoms:
        bonds.update(set(a.bonds))
    for b in bonds:
        a1, a2 = [residue.find_atom(a.name) for a in b.atoms]
        if a1 is None or a2 is None:
            continue
        if a2 not in a1.neighbors:
            m.new_bond(a1, a2)
    from .building import set_new_atom_style
    set_new_atom_style(residue.session, residue.atoms)
    if rename_residue:
        residue.name = template.name

def fix_residue_to_match_md_template(session, residue, md_template, cif_template = None,
        rename=False):
    '''
    For the given residue, add/remove atoms as necessary to match the given MD
    template. If no explicit cif_template argument is given, an attempt will be
    made to find a match using ISOLDE's database. If no CIF template is found,
    corrections to the residue will be limited to deleting excess atoms.
    '''
    import numpy
    from chimerax.isolde.openmm.amberff.template_utils import (
        template_name_to_ccd_name,
        match_template_atoms_to_ccd_atoms
    )
    if cif_template is None:
        ccd_name, _ = template_name_to_ccd_name(md_template.name)
        if ccd_name is not None:
            from chimerax.atomic import mmcif
            cif_template = mmcif.find_template_residue(session, ccd_name)
    else:
        ccd_name = cif_template.name
    if cif_template is None:
        return trim_residue_to_md_template(residue, md_template)
    template_indices, ccd_indices = match_template_atoms_to_ccd_atoms(session, md_template, ccd_name)
    fix_residue_from_template(residue, cif_template, template_indices=ccd_indices)
    template_extra_indices = [i for i in range(len(md_template.atoms)) if i not in template_indices]
    #template_extra_atoms = [md_template.atoms[i] for i in template_extra_indices]
    if len(template_extra_indices):
        template_extra_bonds = set([b for b in md_template.bonds if any([i in template_extra_indices for i in b])])
        from collections import defaultdict
        stub_map = defaultdict(list)
        # stub_map maps an existing atom in the residue to any atoms in the MD
        # template that should be connected to it, but aren't yet modelled.
        found_bonds = set()
        for b in template_extra_bonds:
            i1, i2 = b
            i1_index = numpy.where(template_indices==i1)[0]
            i2_index = numpy.where(template_indices==i2)[0]
            if not len(i1_index) and not len(i2_index):
                continue
            if len(i2_index):
                i1, i2 = i2, i1
                i1_index = i2_index
            i1_index = i1_index[0]
            ccd_atom = cif_template.atoms[ccd_indices[i1_index]]
            res_atom = residue.find_atom(ccd_atom.name)
            if not res_atom:
                raise RuntimeError("Atom {} should be in residue, but isn't".format(ccd_atom.name))
            stub_map[res_atom].append(i2)
            found_bonds.add(b)
        template_extra_bonds = template_extra_bonds.difference(found_bonds)
        from chimerax.core.errors import UserError
        if len(template_extra_bonds):
            err_str = ('MD template {} contains extra atoms that are not in '
                'CCD template {}, and are not directly connected to existing '
                'atoms. Since MD templates do not explicitly provide geometry,'
                'these atoms will not be built. As it stands, the resulting '
                'residue will contain only those atoms which the MD and CCD '
                'templates have in common.').format(md_template.name, ccd_name)
            raise UserError(err_str)
        seen = set()
        for new_atom_list in stub_map.values():
            for i in new_atom_list:
                if i in seen:
                    err_str = ('The atom {} in MD template {} bonds to more than '
                        'one existing atom in residue {}. Since MD templates do '
                        'not explicitly specify geometry, this type of atom addition '
                        'is not currently supported. The resulting residue will '
                        'contain only those atoms which the MD and CCD templates '
                        'have in common').format(
                            md_template.atoms[i].name, md_template.name, ccd_name)
                    raise UserError(err_str)
                seen.add(i)
        from chimerax.atomic import Element
        from chimerax.atomic.build_structure import modify_atom
        for existing_atom, new_indices in stub_map.items():
            num_new_atoms = len(new_indices)
            num_existing_neighbors = len(existing_atom.neighbors)
            num_bonds = len(existing_atom.neighbors) + num_new_atoms
            new_tatoms = [md_template.atoms[i] for i in new_indices]
            from chimerax.atomic.build_structure.mod import ParamError
            try:
                modified_atoms = modify_atom(existing_atom, existing_atom.element,
                    num_bonds, res_name=residue.name)
            except ParamError:
                from chimerax.core.errors import UserError
                err_str = ('Failed to add atoms {} to atom {} because this will '
                    'lead to having {} atoms attached, which is more than its '
                    'assigned geometry can support. This is probably due to an '
                    'error in the MD template ({}). If this template is built '
                    'into ISOLDE, please report this using Help/Report a bug').format(
                    [a.name for a in new_tatoms], existing_atom.name,
                    num_existing_neighbors+len(new_tatoms), md_template.name
                )
                raise UserError(err_str)
            new_atoms = modified_atoms[1:]
            for na, ta in zip(new_atoms, new_tatoms):
                modify_atom(na, Element.get_element(ta.element.atomic_number),
                    1, name=ta.name, res_name=residue.name
                )


def trim_residue_to_md_template(residue, md_template):
    residue.session.logger.status('Trimming residue to MD template...')
    if len(md_template.atoms) > 2:
        from chimerax.isolde.graph import make_graph_from_residue
        rgraph = make_graph_from_residue(residue)
        ri, ti, _ = rgraph.maximum_common_subgraph(md_template.graph, timeout=5)
        if len(ti) != len(md_template.atoms):
            from chimerax.core.errors import UserError
            err_string = ('Template {} contains atoms not found in residue '
            '{} {}{}, and no matching CIF template is available to provide '
            'coordinates.'.format(md_template.name, residue.name, residue.chain_id, residue.number))
            raise UserError(err_string)
        residue.atoms.subtract(residue.atoms[ri]).delete()


def build_next_atom_from_coords(residue, found_neighbors, template_new_atom):
    # print('Building next atom {} from aligned coords of {}'.format(template_new_atom.name,
    #     ','.join([f[0].name for f in found_neighbors])))
    bonded_to = [f[0] for f in found_neighbors]
    m = residue.structure
    n = len(found_neighbors)
    found = set(f[0] for f in found_neighbors)
    while len(found_neighbors) < 3:
        for ra, ta in found_neighbors:
            for tn in ta.neighbors:
                rn = residue.find_atom(tn.name)
                if rn and rn not in found:
                    found_neighbors.append((rn, tn))
                    found.add(rn)
        if len(found_neighbors) <= n:
            residue.session.logger.warning("Couldn't find more than two neighbor atoms. Falling back to simple geometry-based addition.")
            return build_next_atom_from_geometry(residue, *found_neighbors[0], template_new_atom)
        n = len(found_neighbors)
    from chimerax.geometry import align_points
    import numpy
    ra_coords = numpy.array([f[0].coord for f in found_neighbors])
    ta_coords = numpy.array([f[1].coord for f in found_neighbors])
    bfactor = numpy.mean([f[0].bfactor for f in found_neighbors])
    occupancy = numpy.mean([f[0].occupancy for f in found_neighbors])
    tf, rms = align_points(ta_coords, ra_coords)
    from chimerax.atomic.struct_edit import add_atom
    ta = template_new_atom
    a = add_atom(ta.name, ta.element, residue, tf*ta.coord, occupancy=occupancy, bfactor=bfactor)
    for b in bonded_to:
        m.new_bond(a, b)




def build_next_atom_from_geometry(residue, residue_anchor, template_anchor, template_new_atom):
    from chimerax.atomic import struct_edit
    from chimerax.core.geometry import distance, angle, dihedral
    r = residue
    m = r.structure
    tnext = template_new_atom
    if tnext is None:
        raise TypeError('Template does not contain an atom with that name!')
    tstub = template_anchor
    rstub = residue_anchor

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

    # print('Building next atom {} from geometry of {}'.format(template_new_atom.name,
    #     ','.join([n.name for n in (n1, n2, n3)])))

    dist = distance(tnext.coord, tstub.coord)
    ang = angle(tnext.coord, tstub.coord, a2.coord)
    dihe = dihedral(tnext.coord, tstub.coord, a2.coord, a3.coord)
    # print('{}: {} {} {}'.format(next_atom_name, dist, ang, dihe))
    a = struct_edit.add_dihedral_atom(tnext.name, tnext.element, n1, n2, n3, dist, ang, dihe)
    a.occupancy = rstub.occupancy
    a.bfactor = rstub.bfactor
    return a


def copy_ideal_coords_to_exp(ciffile):
    '''
    A temporary measure: ChimeraX currently uses the experimental rather than
    ideal coordinates when loading a template from the CCD. This is problematic
    because the experimental coordinates may be undefined for some atoms. This
    script overwrites the experimental coordinates with the ideal ones for a
    CIF file containing a single residue definition, writing the result as
    {original name}_ideal.cif.
    '''
    import os
    name = os.path.splitext(ciffile)[0]
    with open(ciffile, 'rt') as infile, open(name+'_ideal.cif', 'wt') as outfile:
        lines = infile.read().split('\n')
        count = 0
        line = infile.readline()
        for i, line in enumerate(lines):
            if line.startswith('_chem_comp_atom'):
                break
            if line.startswith('_chem_comp.id'):
                resid = line.split()[-1]
            outfile.write(line+'\n')
        count += i
        exp_indices = [-1,-1,-1]
        ideal_indices=[-1,-1,-1]
        for i, line in enumerate(lines[count:]):
            if not line.startswith('_chem_comp_atom'):
                break
            outfile.write(line+'\n')
            if 'model_Cartn_' in line:
                if not 'ideal' in line:
                    target = exp_indices
                else:
                    target = ideal_indices
                if '_x' in line:
                    target[0] = i
                elif '_y' in line:
                    target[1] = i
                elif '_z' in line:
                    target[2] = i
        # print('Ideal indices: {}, experimental indices: {}'.format(ideal_indices, exp_indices))
        count += i
        for i, line in enumerate(lines[count:]):
            if not line.startswith(resid):
                break
            split_line = line.split()
            for i in range(3):
                split_line[exp_indices[i]] = split_line[ideal_indices[i]]
            outfile.write(' '.join(split_line)+'\n')
        count += i
        for line in lines[count:]:
            outfile.write(line+'\n')

def load_cif_templates(session, cif_files):
    from chimerax import mmcif
    all_names = []
    for cif_file in cif_files:
        try:
            current_names = get_template_names(cif_file, filter_out_obsolete=False)
            if not len(current_names):
                session.logger.warning("File {} does not appear to contain any "
                    "valid residue templates.".format(cif_file))
                session.logger.status("At least one file did not load correctly. Check the log.")
                continue
            mmcif.load_mmCIF_templates(cif_file)
            session.logger.info("Loaded CIF templates for [{}] from {}".format(
                ', '.join(current_names), cif_file
            ))
            all_names.extend(current_names)
        except Exception as e:
            err_string = ("Attempt to load CIF template(s) [{}] from {} failed "
                "with the following error. Are you sure this is a valid template "
                "file? \n {}"
            ).format(', '.join(current_names), cif_file, str(e))
            session.logger.warning(err_string)
            continue
    if not len(all_names):
        raise RuntimeError('No files successfully loaded.')



def get_template_names(cif_file, filter_out_obsolete=True):
    '''
    Get the names of all non-obsolete templates in a chemical components cif file.
    '''
    from chimerax import mmcif
    tables = mmcif.get_cif_tables(cif_file, ['chem_comp','chem_comp_atom'], all_data_blocks=True)
    tables = [t for t in tables if t[1][1].has_field('model_Cartn_x') or t[1][1].has_field('pdbx_model_Cartn_x_ideal')]
    if filter_out_obsolete:
        residue_names = [t[0] for t in tables if t[1][0].fields(['pdbx_release_status'])[0][0] != 'OBS']
    else:
        residue_names = [t[0] for t in tables]
    return residue_names

def match_ff_templates_to_ccd_templates(session, forcefield, ccd_names):
    '''
    For each template in the MD forcefield, try to find the CCD template that
    most closely matches it. This will take a long time!
    '''
    from chimerax.isolde.graph import make_graph_from_residue_template
    from chimerax import mmcif
    ff_templates = forcefield._templates
    best_matches = {}
    ccd_graphs = {}
    small_ccd_templates = []
    for name in ccd_names:
        ccd_template = mmcif.find_template_residue(session, name)
        if len(ccd_template.atoms) > 2:
            ccd_graphs[name] = make_graph_from_residue_template(ccd_template)
        else:
            small_ccd_templates.append(ccd_template)
    for tname, template in ff_templates.items():
        best_score = -1000
        if len(template.atoms) < 4:
            best_match, score = match_small_ff_template_to_ccd_templates(template, small_ccd_templates)
            best_matches[tname] = (best_match, score)
            continue
        tgraph = template.graph
        num_heavy_atoms = len(tgraph.labels[tgraph.labels != 1])
        for ccd_name, ccd_graph in ccd_graphs.items():
            ccd_heavy_atoms = len(ccd_graph.labels[ccd_graph.labels!=1])
            if abs(ccd_heavy_atoms - num_heavy_atoms) > 2:
                continue
            fti, cti, _ = tgraph.maximum_common_subgraph(ccd_graph)
            score = len(fti)*3 - len(tgraph.labels) - len(ccd_graph.labels)
            if score > best_score:
                best_matches[tname] = ([ccd_name], score)
                best_score = score
            elif score == best_score:
                best_matches[tname][0].append(ccd_name)
        print('Best matches for {}: {}'.format(tname, ', '.join(best_matches[tname][0])), flush=True)
    return best_matches

def match_small_ff_template_to_ccd_templates(ff_template, ccd_templates):
    '''
    For a MD template with two or fewer atoms, find any CCD templates with the
    same cohort of elements.
    '''
    found = []
    elements = list(sorted(a.element.atomic_number for a in ff_template.atoms))
    for ct in ccd_templates:
        if len(ct.atoms) != len(ff_template.atoms):
            continue
        c_elements = list(sorted(a.element.number for a in ct.atoms))
        matches = all(cn==n for cn, n in zip(c_elements, elements))
        if matches:
            found.append(ct.name)
    if len(found):
        score = len(elements)
    else:
        score = 0
    return found, score
