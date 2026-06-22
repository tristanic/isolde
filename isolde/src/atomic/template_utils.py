
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

def atom_names_unique(residue):
    return len(residue.atoms) == len(set(residue.atoms.names))

def make_atom_names_unique(residue):
    from collections import defaultdict
    name_map = defaultdict(int)
    for a in residue.atoms:
        name_map[a.element.name] += 1
        a.name = f'{a.element.name}{name_map[a.element.name]:02d}'

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

    This method compares RDKit representations of residue and template to find
    the maximal overlap, returning a dict mapping atoms in the residue to those
    in the template. Matching is bond-order- and chirality-aware: an inverted
    stereocentre (e.g. a D-sugar against an L-template) is excluded from the
    correspondence so it is rebuilt from the template geometry rather than
    silently kept with the wrong coordinates.

    Args:
        * residue
            - a :class:`chimerax.atomic.Residue` object
        * template
            - a :class:`chimerax.atomic.cytmpl.TmplResidue` or a :class:`chimerax.atomic.Residue` object

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
    from chimerax.atomic import Atoms
    if len(residue.atoms) == 1:
        a = residue.atoms[0]
        ta = template.find_atom(a.name)
        if ta is not None:
            return {a: ta}, Atoms(), [a for a in template.atoms if a != ta]
    bonds = residue.atoms.intra_bonds
    if len(bonds) == 0:
        add_bonds_from_template_by_matched_names(residue, template)
    bonds = residue.atoms.intra_bonds
    if len(bonds) == 0:
        raise UserError(f'Residue {residue.name} {residue.chain_id}{residue.number}{residue.insertion_code} has no bonds, and no atom names '
            f'in common with template {template.name}. This will need to be corrected manually.')
    from .template_proxy import TemplateProxy
    from . import rdkit_bridge as rb
    proxy = TemplateProxy(residue.session, template)

    # Convert both sides to RDKit molecules with stereochemistry assigned from
    # coordinates (heavy-atom graphs; H handled separately by name below).
    use_names = (match_by == 'name' and atom_names_unique(residue))
    res_mol, res_map = rb.residue_to_rdkit(residue)
    tmpl_mol, tmpl_map = proxy.to_rdkit()
    # The FMCS path *is* the correspondence engine, so a conversion failure there
    # is fatal (and actionable). The name path only needs RDKit for the chirality
    # check, so it degrades gracefully when conversion fails.
    if not use_names:
        if res_mol is None:
            raise UserError(
                f'Could not interpret the chemistry of residue {residue.name} '
                f'{residue.chain_id}{residue.number}{residue.insertion_code}. It '
                'may be missing atoms or have badly distorted geometry. Try '
                'adding hydrogens, or delete it and rebuild with '
                '"isolde add ligand".')
        if tmpl_mol is None:
            raise UserError(
                f'Could not interpret the chemistry of template {proxy.name}.')

    # Heavy-atom correspondence. In 'name' mode (unique names) trust the names;
    # otherwise use the bond-order-aware maximum common subgraph.
    amap = {}
    if use_names:
        for a in residue.atoms:
            if a.element.number == 1:
                continue
            ta = proxy.find_atom(a.name)
            if ta is not None and ta.element.name == a.element.name:
                amap[a] = ta
    else:
        for ri, ti in rb.fmcs_index_correspondence(res_mol, tmpl_mol):
            ra = res_map.get(ri)
            ta = tmpl_map.get(ti)
            if ra is not None and ta is not None:
                amap[ra] = ta

    # Inverted stereocentres are intentionally kept in the correspondence (they
    # are the *same* atom, just mis-handed). Chirality is repaired separately and
    # locally by correct_chirality_to_template() after the structural fix, which
    # is far less disruptive than deleting and rebuilding a complete centre.

    # Match hydrogens by name (the heavy-atom RDKit graphs carry no H), so any
    # already-correct H aren't deleted as "extra" and missing template H are
    # still flagged for building -- preserving the previous behaviour.
    matched_t = set(amap.values())
    for a in residue.atoms:
        if a.element.number != 1:
            continue
        ta = proxy.find_atom(a.name)
        if ta is not None and ta.element.number == 1 and ta not in matched_t:
            amap[a] = ta
            matched_t.add(ta)

    tatoms = template.atoms
    if limit_template_indices is not None:
        allowed = set(template.atoms[i] for i in limit_template_indices)
        amap = {ra: ta for ra, ta in amap.items() if ta in allowed}
        tatoms = [template.atoms[i] for i in limit_template_indices]

    residue_extra_atoms = Atoms(set(residue.atoms).difference(set(amap.keys())))
    template_extra_atoms = list(set(tatoms).difference(set(amap.values())))

    return amap, residue_extra_atoms, template_extra_atoms

def add_metal_bonds_from_template(residue, template):
    m = residue.structure
    metal_atoms = [a for a in template.atoms if a.element.is_metal]
    from chimerax.atomic.struct_edit import add_bond
    for met in metal_atoms:
        rmet = residue.find_atom(met.name)
        if rmet is None:
            raise TypeError('Residue does not have the required metal ion! Use fix_residue_from_template() instead.')
        for n in met.neighbors:
            rn = residue.find_atom(n.name)
            if not rn in rmet.neighbors:
                add_bond(rmet, rn)


def _place_rigid_blocks(residue, template, missing):
    '''Place as many of the `missing` template atoms as possible by rigid-block
    superposition, returning the template atoms NOT placed (for the caller's
    per-atom fallback).

    The template is decomposed into rigid fragments (ring systems + terminal
    groups; rotatable bonds cut). For each fragment that still has >=3 heavy atoms
    present in the residue, the template fragment is superposed onto those anchors
    (`align_points`) and the fragment's missing atoms -- hydrogens included, riding
    along in the same rigid frame -- are added by that transform. This avoids the
    error accumulation of the per-atom internal-coordinate builder and, because Hs
    travel with their heavy neighbours, places them on the correct side.
    '''
    import numpy
    from . import rdkit_bridge as rb
    session = residue.structure.session
    mol = None
    try:
        mol, _status = rb.template_to_rdkit(session, template.name)
    except Exception:
        mol = None
    if mol is None:
        try:
            mol = rb.residue_to_rdkit(template)
        except Exception:
            mol = None
    if mol is None:
        return set(missing)
    try:
        heavy_groups = rb.rigid_fragments(mol)
    except Exception:
        heavy_groups = []
    if not heavy_groups:
        return set(missing)
    name_to_frag = {}
    for fid, names in enumerate(heavy_groups):
        for nm in names:
            name_to_frag[nm] = fid

    def frag_of(tatom):
        fid = name_to_frag.get(tatom.name)
        if fid is not None:
            return fid
        # Hydrogens (and anything the CCD mol omitted) inherit a heavy neighbour's
        # fragment.
        for nb in tatom.neighbors:
            fid = name_to_frag.get(nb.name)
            if fid is not None:
                return fid
        return None

    from collections import defaultdict
    frag_atoms = defaultdict(list)
    for a in template.atoms:
        fid = frag_of(a)
        if fid is not None:
            frag_atoms[fid].append(a)

    from chimerax.geometry import align_points
    from chimerax.atomic.struct_edit import add_atom
    placed = set()
    for fid, tatoms in frag_atoms.items():
        frag_missing = [a for a in tatoms if a in missing]
        if not frag_missing:
            continue
        # Anchors: this fragment's heavy template atoms already present in the
        # residue (heavy only -- more reliable than possibly-sparse hydrogens).
        anchors = []
        for ta in tatoms:
            if ta.element.number == 1:
                continue
            ra = residue.find_atom(ta.name)
            if ra is not None:
                anchors.append((ta, ra))
        if len(anchors) < 3:
            continue  # leave to the per-atom fallback / 3b torsion placement
        t_coords = numpy.array([ta.coord for ta, ra in anchors])
        r_coords = numpy.array([ra.coord for ta, ra in anchors])
        tf, rms = align_points(t_coords, r_coords)
        bfactor = float(numpy.mean([ra.bfactor for ta, ra in anchors]))
        occupancy = float(numpy.mean([ra.occupancy for ta, ra in anchors]))
        for ta in frag_missing:
            if residue.find_atom(ta.name) is not None:
                continue
            add_atom(
                ta.name,
                ta.element,
                residue,
                tf * ta.coord,
                occupancy=occupancy,
                bfactor=bfactor
            )
            placed.add(ta)
    return set(missing) - placed


def fix_residue_from_template(residue, template, rename_atoms_only=False,
        rename_residue=False, match_by='name', template_indices=None,
        optimise_torsions=True):
    import numpy
    from chimerax.atomic import Atoms
    from chimerax.atomic.struct_edit import add_bond
    if any([numpy.any(numpy.isnan(a.coord)) for a in template.atoms]):
        raise TypeError('Template is missing one or more atom coordinates!')
    # Everything downstream matches and bonds atoms by name. Heavy atoms are safe
    # (matched via the bond-aware graph and renamed to the template's names), but
    # hydrogens are matched purely by name -- and a residue's existing H (e.g. from
    # `addh`) may follow a different convention than the template, so the template's
    # "H21" can sit on a different atom than the residue's "H21". That cross-matches
    # hydrogens, yielding duplicated and cross-bonded H. Strip them up front so
    # every H is (re)built (plain path) or re-added to match the MD template
    # (heavy-only path) with the template's own naming. rename_atoms_only must stay
    # non-destructive, so it keeps its hydrogens.
    if not rename_atoms_only:
        existing_h = residue.atoms[residue.atoms.element_names == 'H']
        if len(existing_h):
            existing_h.delete()
    matched_nodes, residue_extra, template_extra = find_maximal_isomorphous_fragment(
        residue, template, limit_template_indices=template_indices, match_by=match_by
    )

    if len(matched_nodes) < 3:
        from chimerax.atomic import Residue
        if residue.polymer_type == Residue.PT_NONE:
            final_string = (f'Try deleting this residue and replacing it with "isolde add ligand {residue.name} - or have you just forgotten '
            f'to add hydrogens?')
        elif residue.polymer_type == Residue.PT_AMINO:
            final_string = (f'Try deleting this residue and replacing it with "isolde add aa {residue.name}"')
        else:
            final_string = ('Adding of nucleic acid residues is not currently possible in ISOLDE. To move forward, just delete this residue.')

        raise UserError(f'Residue {residue.name} {residue.chain_id}{residue.number}{residue.insertion_code} '
            f'has only {len(matched_nodes)} atoms in common with the template. At least 3 are needed to rebuild automatically. '
            f'{final_string}')

    m = residue.structure
    session = m.session
    rnames = set(residue.atoms.names)
    tnames = set([a.name for a in template.atoms])

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
        return [], []

    # Place whole rigid fragments by superposition first (no error accumulation);
    # anything that can't be block-placed (e.g. a fully-missing fragment) falls
    # through to the per-atom internal-coordinate builder below.
    still_missing = _place_rigid_blocks(residue, template, set(template_extra))
    remaining = set(still_missing)
    while remaining:
        built_this_pass = set()
        for ta in remaining:
            # Candidate anchors: this atom's template neighbours already present in
            # the residue. Try each in turn -- the first whose own connectivity can
            # define a dihedral wins. A neighbour that is itself only just placed
            # (so its further neighbours aren't built yet) is a dead-end anchor and
            # is skipped, not crashed on; the atom is then retried next pass once a
            # better-connected neighbour exists.
            anchors = []
            for tn in ta.neighbors:
                ra = residue.find_atom(tn.name)
                if ra:
                    anchors.append((ra, tn))
            for ra, tn in anchors:
                try:
                    build_next_atom_from_geometry(residue, ra, tn, ta)
                    built_this_pass.add(ta)
                    break
                except InsufficientAnchorsError:
                    continue
        if not built_this_pass:
            raise RuntimeError(
                'Could not rebuild residue {} {}{}: the remaining atoms ({}) have '
                'too few connected reference atoms to place by geometry.'.format(
                    residue.name, residue.chain_id, residue.number,
                    ', '.join(sorted(a.name for a in remaining))
                )
            )
        remaining -= built_this_pass
    bonds = set()
    for a in template.atoms:
        bonds.update(set(a.bonds))
    for b in bonds:
        a1, a2 = [residue.find_atom(a.name) for a in b.atoms]
        if a1 is None or a2 is None:
            continue
        if a2 not in a1.neighbors:
            add_bond(a1, a2)
    # Now that the residue is structurally complete, repair any inverted chiral
    # centres in place (rotate the smallest substituent rather than rebuild).
    corrected, unfixable = [], []
    try:
        corrected, unfixable = correct_chirality_to_template(residue, template)
    except Exception as e:
        session.logger.info(f'Chirality correction skipped for {residue.name}: {e}')
    # Settle freshly-built pendant groups into density / out of clashes by
    # optimising the torsion about each rotatable bond that moves them. The block
    # placement above fixes each fragment's internal geometry but leaves the
    # inter-fragment torsions at whatever the template happened to have.
    if optimise_torsions and len(template_extra):
        try:
            _optimise_pendant_torsions(
                residue, template, set(a.name for a in template_extra))
        except Exception as e:
            session.logger.info(
                f'Pendant-torsion optimisation skipped for {residue.name}: {e}')
    if rename_residue:
        residue.name = template.name
    return corrected, unfixable


def _best_torsion_angle(angles, scores, eps=0.02):
    '''Pick the lowest-scoring angle, but on a (near-)tie prefer the one closest
    to 0 -- i.e. closest to the template/rebuild torsion. Without this, a tie (no
    clashes, no map) would let argmin pick an arbitrary angle and gratuitously
    rotate an already-good pendant. Scores within ``eps`` are treated as tied.
    '''

    def dist_from_zero(a):
        a = a % 360
        return min(a, 360 - a)

    return min(
        zip(angles, scores), key=lambda asc: (round(asc[1] / eps), dist_from_zero(asc[0]))
    )[0]


def _optimise_pendant_torsions(
    residue, template, built_names, map_weight=1.0, clash_weight=1.0
):
    '''Optimise the torsion about each rotatable bond that moves freshly-built
    atoms, to settle pendant groups into density (when an MDFF map is available)
    and out of clashes. A coarse 30-degree scan seeds a local +/-15-degree refine
    per bond -- robust and gradient-free, and a no-op where the template torsion
    is already best (so it never gratuitously disturbs a good rebuild).

    Bonds are taken from the template's rigid-fragment decomposition (an
    inter-fragment, in-residue single bond is rotatable); the moving side is the
    smaller of the two, and only bonds whose moving side contains a freshly-built
    atom are touched. Parent bonds (larger moving side) are optimised before their
    children so a pendant's frame is settled before its own substituents.
    '''
    from . import rdkit_bridge as rb
    from ..refine.fit_utils import near_atoms, pose_score
    model = residue.structure
    session = model.session

    mol = None
    try:
        mol, _status = rb.template_to_rdkit(session, template.name)
    except Exception:
        mol = None
    if mol is None:
        return []
    try:
        heavy_groups = rb.rigid_fragments(mol)
    except Exception:
        return []
    if len(heavy_groups) < 2:
        return []
    name_to_frag = {}
    for fid, names in enumerate(heavy_groups):
        for nm in names:
            name_to_frag[nm] = fid

    # Collect rotatable bonds: in-residue single bonds between heavy atoms in two
    # different rigid fragments, whose smaller side carries a freshly-built atom.
    jobs = []
    for b in residue.atoms.intra_bonds:
        a1, a2 = b.atoms
        if a1.element.number == 1 or a2.element.number == 1:
            continue
        f1 = name_to_frag.get(a1.name)
        f2 = name_to_frag.get(a2.name)
        if f1 is None or f2 is None or f1 == f2:
            continue
        try:
            side1 = b.side_atoms(a1)
            side2 = b.side_atoms(a2)
        except Exception:
            # In a ring/cycle -- not actually rotatable; skip.
            continue
        # Move the side carrying the freshly-built atoms, leaving the kept,
        # correctly-placed side alone -- even if the built side is the larger one
        # (that is exactly the lever that swings a whole diverged assembly back
        # into density). Only when both sides are freshly built does "move the
        # smaller side" apply, to keep the swing local. A bond with no built atom
        # on either side is already settled and is skipped.
        side1_built = any(a.name in built_names for a in side1)
        side2_built = any(a.name in built_names for a in side2)
        if side1_built and side2_built:
            moving_atom, moving_side = (a1, side1) if len(side1) <= len(side2) \
                else (a2, side2)
        elif side1_built:
            moving_atom, moving_side = a1, side1
        elif side2_built:
            moving_atom, moving_side = a2, side2
        else:
            continue
        jobs.append((b, moving_atom, len(moving_side)))

    if not jobs:
        return []
    # Parent-first: larger moving side before smaller.
    jobs.sort(key=lambda j: j[2], reverse=True)

    from chimerax.atomic import Atoms
    from chimerax.geometry import rotation
    from ..refine.fit_utils import active_mdff_map
    volume, _gk = active_mdff_map(model)

    optimised = []
    coarse = list(range(0, 360, 30))
    for b, moving_atom, _size in jobs:
        moving_side = b.side_atoms(moving_atom)
        # The moving-side atom sits on the rotation axis, so axis/centre are fixed
        # across all trial angles; snapshot once and apply each angle absolutely.
        centre = moving_atom.coord
        axis = centre - b.other_atom(moving_atom).coord
        orig = moving_side.coords.copy()
        # Clash partners must cover the WHOLE rotational sweep, not just the start
        # pose -- otherwise a pose that swings the group into nearby protein isn't
        # penalised (the offending atoms were never in the context set). Union the
        # neighbours found at every coarse angle.
        near = set()
        for ang in coarse:
            moving_side.coords = rotation(axis, ang, centre).transform_points(orig)
            near.update(near_atoms(moving_side, model, 5.0))
        moving_side.coords = orig
        surrounds = Atoms(list(near))

        def score_at(angle):
            tf = rotation(axis, angle, centre)
            moving_side.coords = tf.transform_points(orig)
            return pose_score(
                session,
                moving_side,
                surrounds,
                volume,
                map_weight=map_weight,
                clash_weight=clash_weight
            )

        best = _best_torsion_angle(coarse, [score_at(a) for a in coarse])
        fine = [best + d for d in range(-15, 16, 5)]
        best = _best_torsion_angle(fine, [score_at(a) for a in fine])
        # Leave the moving side at the best pose.
        score_at(best)
        if best % 360 != 0:
            optimised.append((moving_atom.name, best % 360))
    return optimised


# Largest substituent (heavy-atom count) we'll rigidly rotate to invert a chiral
# centre. A methyl is 1; small groups (hydroxyl, amino, halide) are fine. Bigger
# rigid groups swing too far and are left for delete-and-rebuild instead. Tunable.
MAX_ROTATABLE_CHIRAL_BRANCH = 4


def _branch_atoms(chiral_atom, root):
    '''Atoms of the substituent branch rooted at `root` (a neighbour of
    `chiral_atom`), traversing within the residue but never back through
    `chiral_atom`.

    Returns ``(atoms, cyclic)`` where ``cyclic`` is True if the branch loops back
    to another neighbour of the chiral atom -- i.e. it is part of a ring and so
    cannot be rigidly rotated to invert the centre.
    '''
    residue = chiral_atom.residue
    # A branch is cyclic if it reaches ANOTHER neighbour of the chiral atom (a
    # ring back through the centre). Reaching the root again is NOT cyclic -- it
    # is just the branch's own internal bonding (e.g. a methyl's hydrogens are
    # bonded to the root carbon), so the root must be excluded here.
    ring_closers = set(chiral_atom.neighbors)
    ring_closers.discard(root)
    seen = {chiral_atom, root}
    out = [root]
    stack = [root]
    cyclic = False
    while stack:
        a = stack.pop()
        for nb in a.neighbors:
            if nb is chiral_atom:
                continue
            if nb in ring_closers:
                cyclic = True       # branch reconnects to the centre via a ring
            if nb in seen or nb.residue is not residue:
                continue
            seen.add(nb)
            out.append(nb)
            stack.append(nb)
    return out, cyclic


def invert_chiral_center(chiral_atom):
    '''Invert the stereochemistry at `chiral_atom` in place by rigidly rotating
    its smallest exocyclic substituent branch 180 degrees about an axis through
    the atom.

    This flips the centre's handedness without *mirroring* the moved branch, so
    the branch's internal geometry is preserved (a reflection could introduce new
    stereo errors in a substituent that is itself chiral). The result is
    geometrically strained by design; energy minimisation relaxes it into the
    correct (now inverted) well, which it cannot escape without crossing the
    planar inversion barrier.

    Returns True if the centre was inverted, or False if it has no rotatable
    (exocyclic) substituent -- e.g. an all-ring centre -- in which case the
    caller should fall back to delete-and-rebuild.
    '''
    import numpy
    from chimerax.geometry import rotation
    neighbors = list(chiral_atom.neighbors)
    if len(neighbors) < 3:
        return False
    # Rank rotatable (acyclic) substituent branches by heavy-atom count, then mass
    # (so a lone hydrogen is preferred, then a methyl, etc.).
    candidates = []
    for nb in neighbors:
        atoms, cyclic = _branch_atoms(chiral_atom, nb)
        if cyclic:
            continue
        n_heavy = sum(1 for a in atoms if a.element.number != 1)
        weight = sum(a.element.mass for a in atoms)
        candidates.append((n_heavy, weight, nb, atoms))
    if not candidates:
        return False
    candidates.sort(key=lambda c: (c[0], c[1]))
    n_heavy, _, root, branch_atoms = candidates[0]
    if n_heavy > MAX_ROTATABLE_CHIRAL_BRANCH:
        # Even the smallest movable substituent is too large to rotate safely;
        # leave this centre for delete-and-rebuild.
        return False

    C = numpy.array(chiral_atom.coord)
    v = numpy.array(root.coord) - C
    # A 180 deg rotation about any axis through C perpendicular to the C-root
    # bond sends the root to the antipodal side, inverting the centre. Bias the
    # axis away from the other substituents where possible to reduce clashes.
    others = [numpy.array(nb.coord) - C for nb in neighbors if nb is not root]
    ref = numpy.sum(others, axis=0) if others else numpy.array([1.0, 0.0, 0.0])
    axis = numpy.cross(v, ref)
    if numpy.linalg.norm(axis) < 1e-6:
        for trial in ([1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]):
            axis = numpy.cross(v, trial)
            if numpy.linalg.norm(axis) >= 1e-6:
                break
    axis = axis / numpy.linalg.norm(axis)
    tf = rotation(axis, 180.0, center=chiral_atom.coord)
    for a in branch_atoms:
        a.coord = tf * a.coord
    return True


def correct_chirality_to_template(residue, template):
    '''Repair inverted chiral centres in `residue` so their handedness matches
    `template`, in place.

    Each mis-handed centre is fixed locally by :func:`invert_chiral_center`
    (rotating its smallest substituent) rather than by deleting and rebuilding --
    appropriate when the residue is structurally complete and only the
    configuration is wrong. Relaxation of the (intentionally strained) result is
    left to the caller's normal energy minimisation / simulation.

    Returns ``(corrected, unfixable)``: lists of the chiral-atom names that were
    corrected in place and that could not be corrected (e.g. all-ring centres with
    no rotatable substituent). Both empty if nothing needed fixing.
    '''
    import numpy
    from . import rdkit_bridge as rb
    from .template_proxy import TemplateProxy
    session = residue.structure.session
    proxy = TemplateProxy(session, template)
    tmpl_mol, _ = proxy.to_rdkit()
    if tmpl_mol is None:
        return [], []
    # The clean template tells us WHICH atoms are stereocentres (reliable CIP).
    tmpl_cip = rb.cip_codes_by_name(tmpl_mol)
    if not tmpl_cip:
        return [], []

    def _signed_vol(center, subs):
        c = numpy.asarray(center)
        v = [numpy.asarray(s) - c for s in subs]
        return float(numpy.dot(v[0], numpy.cross(v[1], v[2])))

    corrected = []
    unfixable = []
    for name in tmpl_cip:
        atom = residue.find_atom(name)
        tatom = template.find_atom(name)
        if atom is None or tatom is None:
            continue
        # Compare handedness by GEOMETRY -- the signed chiral volume of three
        # common heavy substituents, in the SAME order, in the residue vs the
        # template. This is robust to stripped hydrogens and to mis-perceived
        # bonds (e.g. a malformed linkage typed aromatic), unlike comparing RDKit
        # CIP codes of the residue, which can flip a correctly-placed centre when
        # the residue is mid-rebuild or chemically mis-perceived.
        pairs = []
        for tn in tatom.neighbors:
            if tn.element.number == 1:
                continue
            ra = residue.find_atom(tn.name)
            if ra is not None and ra.element.number > 1:
                pairs.append((ra, tn))
            if len(pairs) == 3:
                break
        if len(pairs) < 3:
            continue
        rv = _signed_vol(atom.coord, [ra.coord for ra, tn in pairs])
        tv = _signed_vol(tatom.coord, [tn.coord for ra, tn in pairs])
        if rv == 0.0 or tv == 0.0:
            continue
        if (rv > 0) != (tv > 0):
            if invert_chiral_center(atom):
                corrected.append(name)
            else:
                unfixable.append(name)
    return corrected, unfixable


def add_bonds_from_template_by_matched_names(residue, template):
    from chimerax.atomic.struct_edit import add_bond
    amap = {}
    for a in residue.atoms:
        ta = template.find_atom(a.name)
        if ta is not None:
            amap[a] = ta
    for a, ta in amap.items():
        for n in ta.neighbors:
            rn = residue.find_atom(n.name)
            if rn is not None and rn not in a.neighbors:
                add_bond(a, rn)


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
            from chimerax import mmcif
            try:
                cif_template = mmcif.find_template_residue(session, ccd_name)
            except ValueError:
                pass
    else:
        ccd_name = cif_template.name
    if cif_template is None:
        return trim_residue_to_md_template(residue, md_template)
    template_indices, ccd_indices = match_template_atoms_to_ccd_atoms(session, md_template, ccd_name)
    fix_residue_from_template(residue, cif_template, template_indices=ccd_indices)
    from chimerax.atomic import Atoms
    ratoms = Atoms([residue.find_atom(cif_template.atoms[i].name) for i in ccd_indices])
    residue_indices = residue.atoms.indices(ratoms)
    #template_extra_atoms = [md_template.atoms[i] for i in template_extra_indices]
    add_missing_md_template_atoms(session, residue, md_template, residue_indices, template_indices)

def add_missing_md_template_atoms(session, residue, md_template, residue_indices, template_indices):
    import numpy
    template_extra_indices = [i for i in range(len(md_template.atoms)) if i not in template_indices]
    if not len(template_extra_indices):
        return
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
        res_atom = residue.atoms[residue_indices[i1_index]]
        # if not res_atom:
        #     raise RuntimeError("Atom {} should be in residue, but isn't".format(ccd_atom.name))
        stub_map[res_atom].append(i2)
        found_bonds.add(b)
    template_extra_bonds = template_extra_bonds.difference(found_bonds)
    if len(template_extra_bonds):
        err_str = ('MD template {} for residue {} {}{}{} contains extra atoms that are not in '
            'a coordinate template, and are not directly connected to existing '
            'atoms. Since MD templates do not explicitly provide geometry,'
            'these atoms will not be built.').format(md_template.name,
                residue.name, residue.chain_id, residue.number, residue.insertion_code)
        session.logger.warning(err_str)
    seen = set()
    for new_atom_list in stub_map.values():
        for i in new_atom_list:
            if i in seen:
                err_str = ('The atom {} in MD template {} bonds to more than '
                    'one existing atom in residue {}. Since MD templates do '
                    'not explicitly specify geometry, this type of atom addition '
                    'is not currently supported. The resulting residue will '
                    'contain only those atoms which the MD and coordinate templates '
                    'have in common').format(
                        md_template.atoms[i].name, md_template.name, residue.name)
                raise UserError(err_str)
            seen.add(i)
    from chimerax.atomic import Element
    from chimerax.build_structure import modify_atom
    for existing_atom, new_indices in stub_map.items():
        num_new_atoms = len(new_indices)
        num_existing_neighbors = len(existing_atom.neighbors)
        num_bonds = len(existing_atom.neighbors) + num_new_atoms
        new_tatoms = [md_template.atoms[i] for i in new_indices]
        from chimerax.build_structure.mod import ParamError
        try:
            modified_atoms = modify_atom(existing_atom, existing_atom.element,
                num_bonds, res_name=residue.name)
        except ParamError:
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
    err_string = ('Auto-correction to a MD template is currently only '
            'possible when both residue and template have at least 3 connected atoms in common. '
            'For smaller residues it is best to just delete and replace the '
            'existing one. If your residue is larger than this, double-check to '
            'see if you have chosen the correct MD template.')
    residue.session.logger.status('Trimming residue to MD template...')
    if len(md_template.atoms) > 2 and len(md_template.bonds) > 1 and len(residue.atoms) > 2 and len(residue.atoms.intra_bonds) > 1:
        from chimerax.isolde.graph import make_graph_from_residue
        rgraph = make_graph_from_residue(residue)
        ri, ti, _ = rgraph.maximum_common_subgraph(md_template.graph, big_first=True, timeout=5)
        # if len(ti) != len(md_template.atoms):
        #     from chimerax.core.errors import UserError
        #     err_string = ('Template {} contains atoms not found in residue '
        #     '{} {}{}, and no matching CIF template is available to provide '
        #     'coordinates.'.format(md_template.name, residue.name, residue.chain_id, residue.number))
        #     raise UserError(err_string)
        ratoms = residue.atoms[ri]
        if len(ratoms) < 3 or len(ratoms.intra_bonds) < 2:
            raise UserError(err_string)

        residue.atoms.subtract(ratoms).delete()
        # Need to redo the graph matching, because indices may have changed
        rgraph = make_graph_from_residue(residue)
        ri, ti, _ = rgraph.maximum_common_subgraph(md_template.graph, big_first=True, timeout=5)
        add_missing_md_template_atoms(residue.session, residue, md_template, ri, ti)
    else:
        raise UserError(err_string)



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
    from chimerax.atomic.struct_edit import add_bond
    for b in bonded_to:
        add_bond(a, b)

_geometry_to_angle = {
    3:  180,
    4:  120,
}


class InsufficientAnchorsError(Exception):
    '''Raised by :func:`build_next_atom_from_geometry` when the chosen anchor does
    not yet have enough connected, already-present heavy atoms to define the
    dihedral needed to place the new atom. The caller should try another anchor or
    defer the atom to a later pass -- it is not a hard failure.'''
    pass


def build_next_atom_from_geometry(residue, residue_anchor, template_anchor, template_new_atom):
    from chimerax.atomic import struct_edit
    from chimerax.geometry import distance, angle, dihedral
    r = residue
    m = r.structure
    tnext = template_new_atom
    if tnext is None:
        raise TypeError('Template does not contain an atom with that name!')
    tstub = template_anchor
    rstub = residue_anchor
    existing_rstub_neighbors = rstub.neighbors

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
    elif len(t_direct_neighbors) == 1:
        a2 = t_direct_neighbors[0]
        n2 = r_direct_neighbors[0]
    else:
        raise InsufficientAnchorsError(
            'No n2 found - Not enough connected atoms to form a dihedral!'
        )
    if not n2:
        raise InsufficientAnchorsError(
            'No n2 found - Not enough connected atoms to form a dihedral!'
        )
    if not n3:
        for a3 in a2.neighbors:
            if a3 not in (a2, tstub) and a3.element.name != 'H':
                n3 = r.find_atom(a3.name)
                if n3:
                    break
        if not n3:
            raise InsufficientAnchorsError(
                'No n3 found - Not enough connected atoms to form a dihedral!'
            )

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
