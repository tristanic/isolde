# @Author: Tristan Croll
# @Date:   27-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 27-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Agent-specific REST helpers (the transport/introspection plane).

These are *not* ISOLDE commands and would be noise as interactive CLI verbs:
they exist only for a programmatic/remote client (the MCP server). Per the
command-vs-helper rubric, ID-returning loaders live in :mod:`..server_methods`;
this module holds the agent transport surface:

  - ``agent_tools``   : the auto-generated tool manifest (agent-safe commands)
  - ``agent_invoke``  : invoke one agent-safe command from a JSON args object
  - ``list_models``   : inventory of open models so an eyeless agent gets IDs
  - ``resolve_spec``  : describe what an atom-spec selects, to confirm before acting

Each function takes the ChimeraX session first and returns a JSON-serializable
dict (the REST layer wraps them with thread_safe so they run on the UI thread).
'''


def agent_tools(session):
    '''
    Return the manifest of agent-reachable ISOLDE commands: one typed tool
    definition per command flagged ``agent_safe`` in the agent manifest.
    '''
    from chimerax.isolde.cmd.command_registry import agent_safe_commands
    from chimerax.isolde.cmd.schema import command_tool_definition
    tools = []
    for name, record in sorted(agent_safe_commands().items()):
        try:
            tools.append(command_tool_definition(record))
        except Exception as e:
            tools.append({'name': name.replace(' ', '_'), 'command': name,
                          'error': 'schema generation failed: %s' % e})
    from chimerax.isolde.cmd.agent_philosophy import MANIFEST_PREAMBLE
    return {'preamble': MANIFEST_PREAMBLE, 'tools': tools, 'count': len(tools)}


def agent_philosophy(session):
    '''Return ISOLDE's design philosophy for agents (the single source of truth in
    chimerax.isolde.cmd.agent_philosophy): the condensed ``instructions`` string
    for the MCP server's initialize handshake, the manifest ``preamble``, the full
    ``principles`` list, and a rendered ``document``. The MCP server fetches this
    once at start-up.'''
    from chimerax.isolde.cmd.agent_philosophy import as_payload
    return as_payload()


def agent_invoke(session, name, args=None):
    '''
    Invoke the agent-safe ISOLDE command *name* with a JSON *args* object.

    Returns ``{command, result, witness}`` -- where ``witness`` is the
    post-condition movement/heartbeat proof for mutating commands (or null).
    Raises if the command is unknown or not agent_safe.

    Args:
        name: full command name, e.g. 'isolde sim' or its underscore form 'isolde_sim'
        args: dict of arguments validated against the command's input schema
    '''
    from chimerax.isolde.cmd.invoke import invoke_command
    if args is None:
        args = {}
    # Accept both 'isolde sim' and the MCP underscore form 'isolde_sim'.
    name = _canonical_command_name(name)
    return invoke_command(session, name, args)


def _canonical_command_name(name):
    from chimerax.isolde.cmd.command_registry import get_command
    if get_command(name) is not None:
        return name
    spaced = name.replace('_', ' ')
    if get_command(spaced) is not None:
        return spaced
    return name  # let invoke_command raise a clear KeyError


def list_models(session):
    '''
    Inventory of open models, so an agent (which cannot see the Models panel)
    can discover what is loaded and obtain the ID strings it must reference.
    '''
    from chimerax.atomic import AtomicStructure
    try:
        from chimerax.map import Volume
    except Exception:
        Volume = ()
    isolde = getattr(session, 'isolde', None)
    selected = getattr(isolde, 'selected_model', None) if isolde else None
    selected_id = selected.id_string if selected is not None else None

    models = []
    for m in session.models.list():
        entry = {
            'id': m.id_string,
            'spec': getattr(m, 'atomspec', '#' + m.id_string),
            'name': m.name,
            'type': type(m).__name__,
        }
        if isinstance(m, AtomicStructure):
            entry['kind'] = 'atomic'
            entry['n_atoms'] = m.num_atoms
            try:
                entry['chains'] = sorted(set(m.residues.chain_ids.tolist()))
            except Exception:
                pass
            entry['is_selected'] = (m.id_string == selected_id)
        elif Volume and isinstance(m, Volume):
            entry['kind'] = 'map'
        models.append(entry)
    return {'models': models, 'count': len(models),
            'isolde_selected_model': selected_id}


def _structure_from_spec(session, model):
    '''Resolve a model id/spec to a single AtomicStructure.'''
    from chimerax.atomic import AtomicStructure
    from chimerax.core.commands import AtomSpecArg
    parsed, _t, _r = AtomSpecArg.parse(str(model), session)
    structs = [m for m in parsed.evaluate(session).models
               if isinstance(m, AtomicStructure)]
    if not structs:
        raise ValueError('no atomic structure matches %r' % model)
    return structs[0]


def describe_model(session, model, max_sequence=4000):
    '''
    Structured description of an atomic model for an agent that cannot see it:
    chains with one-letter sequences and modelled/unmodelled state, chain breaks,
    ligands/non-standard residues, and water count.

    Args:
        model: model id or atom-spec (e.g. '#1.2' or '1.2'), as returned by list_models.
        max_sequence: cap per-chain sequence length in the response (truncation flagged).
    '''
    from chimerax.atomic import Residue
    s = _structure_from_spec(session, model)

    pt_name = {getattr(Residue, 'PT_AMINO', 1): 'protein',
               getattr(Residue, 'PT_NUCLEIC', 2): 'nucleic'}

    chains = []
    for ch in s.chains:
        residues = list(ch.residues)            # Residue or None, aligned to characters
        modelled = [r is not None for r in residues]
        n_model = sum(modelled)
        # Unmodelled / break regions: runs of None (gaps in modelled structure),
        # plus numbering discontinuities between adjacent existing residues.
        breaks = []
        run_start = None
        for i, present in enumerate(modelled):
            if not present and run_start is None:
                run_start = i
            elif present and run_start is not None:
                breaks.append(_gap(residues, run_start, i, 'unmodelled'))
                run_start = None
        if run_start is not None:
            breaks.append(_gap(residues, run_start, len(residues), 'unmodelled'))
        # numbering jumps among consecutive existing residues (no None between)
        existing = ch.existing_residues
        for a, b in zip(existing[:-1], existing[1:]):
            if b.number - a.number > 1 and modelled and _no_none_between(residues, a, b):
                breaks.append({'kind': 'numbering',
                               'prev_modelled_resnum': a.number,
                               'next_modelled_resnum': b.number,
                               'gap_length': b.number - a.number - 1})

        seq = ch.characters
        truncated = len(seq) > max_sequence
        chains.append({
            'chain_id': ch.chain_id,
            'polymer_type': pt_name.get(int(ch.polymer_type), 'other'),
            'full_sequence_known': bool(ch.full_sequence_known),
            'num_residues': int(ch.num_residues),
            'num_modelled': int(n_model),
            'num_unmodelled': int(ch.num_residues - n_model),
            'sequence': seq[:max_sequence],
            'sequence_truncated': truncated,
            'chain_breaks': breaks,
        })

    # Ligands / non-standard residues and waters (non-polymer residues).
    waters = 0
    ligands = []
    ligand_counts = {}
    for r in s.residues:
        if int(r.polymer_type) != int(getattr(Residue, 'PT_NONE', 0)):
            continue
        if r.name in ('HOH', 'WAT', 'DOD'):
            waters += 1
            continue
        ligand_counts[r.name] = ligand_counts.get(r.name, 0) + 1
        if len(ligands) < 200:
            ligands.append({'resname': r.name, 'chain_id': r.chain_id,
                            'number': r.number, 'spec': r.string(style='command'),
                            'num_atoms': r.num_atoms})

    return {
        'model': s.atomspec,
        'id': s.id_string,
        'name': s.name,
        'num_atoms': s.num_atoms,
        'num_residues': s.num_residues,
        'num_chains': len(chains),
        'chains': chains,
        'ligands': ligands,
        'ligand_summary': ligand_counts,
        'num_waters': waters,
    }


def _gap(residues, start, end, kind):
    '''Describe a run residues[start:end] of unmodelled positions by the modelled
    residue numbers flanking it. A None side means the gap is terminal (no
    modelled residue on that side, e.g. an unmodelled N- or C-terminus).'''
    prev_n = residues[start - 1].number if start > 0 and residues[start - 1] else None
    next_n = residues[end].number if end < len(residues) and residues[end] else None
    return {'kind': kind, 'prev_modelled_resnum': prev_n,
            'next_modelled_resnum': next_n, 'gap_length': end - start}


def _no_none_between(residues, a, b):
    try:
        ia, ib = residues.index(a), residues.index(b)
        return not any(residues[i] is None for i in range(ia + 1, ib))
    except ValueError:
        return True


def resolve_spec(session, spec):
    '''
    Resolve a ChimeraX atom-spec string and describe what it selects (counts +
    bounding box + centroid), so an agent can confirm a selection is what it
    intended before issuing a mutating command against it.
    '''
    from chimerax.core.commands import AtomSpecArg
    try:
        parsed, _text, _rest = AtomSpecArg.parse(spec, session)
        objects = parsed.evaluate(session)
    except Exception as e:
        return {'spec': spec, 'valid': False, 'error': str(e)}

    atoms = objects.atoms
    out = {
        'spec': spec,
        'valid': True,
        'n_atoms': len(atoms),
        'n_residues': objects.num_residues if hasattr(objects, 'num_residues') else None,
        'n_models': len(objects.models),
        'model_ids': [m.id_string for m in objects.models],
    }
    if len(atoms):
        try:
            import numpy
            coords = atoms.scene_coords
            out['centroid'] = coords.mean(axis=0).tolist()
            out['bbox_min'] = coords.min(axis=0).tolist()
            out['bbox_max'] = coords.max(axis=0).tolist()
        except Exception:
            pass
        # A few representative residue specs for orientation.
        try:
            res = atoms.unique_residues
            out['residue_specs'] = [r.string(style='command') for r in res[:25]]
            out['n_unique_residues'] = len(res)
        except Exception:
            pass
    return out


# ---------------------------------------------------------------------------
# Perception helpers (Phase F): where-to-focus, restraint state, per-residue
# detail, density fit, B-factor outliers, view/selection state.
# ---------------------------------------------------------------------------

def problem_zones(session, model, outliers_only=True, cutoff=3, min_points=6):
    '''
    Cluster the current geometry/restraint problems in a model into spatial
    "zones" so the agent knows where to focus. Aggregates unsatisfied restraints
    + validation outliers (rotamers, backbone/Rama, clashes) and DBSCAN-clusters
    them. Returns clustered zones (with centroid + member residues + a type
    breakdown) and an unclustered remainder.
    '''
    s = _structure_from_spec(session, model)
    from chimerax.isolde.problem_regions.problems import ProblemAggregator
    import numpy
    agg = ProblemAggregator(session)
    clustered, remainder = agg.problem_zones(
        s, cutoff=cutoff, min_points=min_points,
        validation_outliers_only=outliers_only)

    def _norm(items):
        # itemgetter(*idx) returns a bare item when a cluster has length 1.
        return list(items) if isinstance(items, (list, tuple)) else [items]

    def _types(items):
        from collections import Counter
        c = Counter()
        for site in items:
            try:
                c[agg.registered_name(type(site))] += 1
            except Exception:
                c[type(site).__name__] += 1
        return dict(c)

    def _residues(items):
        try:
            return [r.string(style='command')
                    for r in agg.cluster_atoms(items).unique_residues[:30]]
        except Exception:
            return []

    zones = []
    for raw in clustered:
        items = _norm(raw)
        try:
            center = numpy.array([site.center for site in items]).mean(axis=0).tolist()
        except Exception:
            center = None
        zones.append({'size': len(items), 'center': center,
                      'types': _types(items), 'residues': _residues(items)})
    zones.sort(key=lambda z: z['size'], reverse=True)
    rem = _norm(remainder) if len(remainder) else []
    return {'model': s.atomspec, 'n_zones': len(zones), 'zones': zones,
            'n_unclustered': len(rem), 'unclustered_types': _types(rem)}


def restraint_summary(session, model):
    '''Per restraint-type counts (total / enabled / unsatisfied) for a model, so
    the agent can see what restraints are active and whether they are satisfied.'''
    s = _structure_from_spec(session, model)
    from chimerax.isolde import session_extensions as sx
    out = {}

    residues = getattr(s, 'residues', None)

    def _collection(mgr):
        coll = getattr(mgr, 'all_restraints', None)
        if coll is not None:
            return coll
        getter = getattr(mgr, 'get_all_restraints_for_residues', None)
        if getter is not None and residues is not None:
            try:
                return getter(residues)
            except Exception:
                return None
        return None

    def _tally(label, mgr):
        if mgr is None:
            return
        coll = _collection(mgr)
        if coll is None:
            try:
                out[label] = {'total': int(mgr.num_restraints)}
            except Exception:
                pass
            return
        try:
            enabled = coll[coll.enableds]
            n_enabled = len(enabled)
            n_unsat = int(len(enabled[enabled.unsatisfieds])) if n_enabled else 0
            out[label] = {'total': len(coll), 'enabled': n_enabled,
                          'unsatisfied': n_unsat}
        except Exception:
            out[label] = {'total': len(coll)}

    for getter, label in (
            ('get_distance_restraint_mgr', 'distance'),
            ('get_proper_dihedral_restraint_mgr', 'proper_dihedral'),
            ('get_adaptive_dihedral_restraint_mgr', 'adaptive_dihedral'),
            ('get_position_restraint_mgr', 'position'),
            ('get_chiral_restraint_mgr', 'chiral')):
        try:
            _tally(label, getattr(sx, getter)(s, create=False))
        except Exception:
            pass
    try:
        for mgr in sx.get_all_adaptive_distance_restraint_mgrs(s):
            _tally('adaptive_distance:' + getattr(mgr, 'name', '?'), mgr)
    except Exception:
        pass
    return {'model': s.atomspec, 'restraints': out}


def _ss_name(residue):
    try:
        if residue.is_helix:
            return 'helix'
        if residue.is_strand:
            return 'strand'
    except Exception:
        return None
    return 'coil'


def residue_info(session, spec, limit=50):
    '''Per-residue detail for targeted fixes: rotamer score/options, Ramachandran
    case/score/phi-psi, secondary structure, and B-factor stats.'''
    import numpy
    from chimerax.core.commands import AtomSpecArg
    from chimerax.isolde.session_extensions import (
        get_ramachandran_mgr, get_rotamer_mgr)
    parsed, _t, _r = AtomSpecArg.parse(str(spec), session)
    residues = parsed.evaluate(session).atoms.unique_residues
    rama_mgr = get_ramachandran_mgr(session)
    rota_mgr = get_rotamer_mgr(session)
    try:
        rama_by_res = {ra.residue: ra for ra in rama_mgr.get_ramas(residues)}
    except Exception:
        rama_by_res = {}
    try:
        rota_by_res = {ro.residue: ro for ro in rota_mgr.get_rotamers(residues)}
    except Exception:
        rota_by_res = {}

    def deg(rad):
        return float(numpy.degrees(rad))

    items = []
    for r in residues[:limit]:
        info = {'spec': r.string(style='command'), 'name': r.name,
                'number': r.number, 'chain': r.chain_id,
                'secondary_structure': _ss_name(r)}
        try:
            b = r.atoms.bfactors
            info['bfactor'] = {'mean': float(b.mean()), 'min': float(b.min()),
                               'max': float(b.max())}
        except Exception:
            pass
        ra = rama_by_res.get(r)
        if ra is not None:
            try:
                if ra.valid:
                    pp = ra.phipsi
                    info['rama'] = {'case': int(ra.case), 'score': float(ra.score),
                                    'phi': deg(pp[0]), 'psi': deg(pp[1])}
            except Exception:
                pass
        ro = rota_by_res.get(r)
        if ro is not None:
            try:
                nt = ro.nearest_target
                info['rotamer'] = {'score': float(ro.score),
                                   'num_targets': int(ro.num_targets),
                                   'nearest_target': nt.get('name') if isinstance(nt, dict) else None}
            except Exception:
                pass
        items.append(info)
    return {'residues': items, 'count': len(items),
            'truncated': len(residues) > limit}


def map_info(session, model):
    '''Density / map-fit query. For crystallographic models reports Rwork/Rfree
    (computed by Clipper each map recalc); for cryo-EM there is no reliable scalar
    metric. Also reports MDFF coupling state. See fit_guidance.'''
    s = _structure_from_spec(session, model)
    try:
        from chimerax.clipper.symmetry import get_map_mgr
    except Exception:
        return {'model': s.atomspec, 'has_maps': False,
                'note': 'chimerax.clipper not available'}
    mmgr = get_map_mgr(s, create=False)
    if mmgr is None:
        return {'model': s.atomspec, 'has_maps': False,
                'note': 'no maps associated with this model'}

    out = {'model': s.atomspec, 'has_maps': True,
           'crystallographic': [], 'cryo_em': None, 'mdff': []}
    for xset in mmgr.xmapsets:
        entry = {'name': getattr(xset, 'name', None), 'rwork': None,
                 'rfree': None, 'resolution': None}
        for attr in ('rwork', 'rfree', 'resolution'):
            try:
                entry[attr] = float(getattr(xset, attr))
            except Exception:
                pass
        try:
            entry['maps'] = [m.name for m in xset.all_maps]
        except Exception:
            pass
        out['crystallographic'].append(entry)
    try:
        nx = mmgr.nxmapset
        if nx is not None:
            out['cryo_em'] = {'maps': [m.name for m in nx.all_maps]}
    except Exception:
        pass

    # MDFF coupling state per map (read-only here; see set_mdff to change).
    from chimerax.isolde.session_extensions import get_mdff_mgr
    try:
        for v in mmgr.all_maps:
            mgr = get_mdff_mgr(s, v, create=False)
            if mgr is not None:
                out['mdff'].append({'volume': v.name,
                                    'enabled': bool(mgr.enabled),
                                    'coupling_constant': float(mgr.global_k)})
    except Exception:
        pass

    out['fit_guidance'] = (
        'Crystallographic: watch Rfree/Rwork, but they can rise slightly after a '
        'genuine improvement pending B-factor re-refinement; the difference map is '
        'also informative. Cryo-EM: no reliable scalar fit metric exists -- use '
        'render for visual confirmation. Do NOT raise MDFF coupling to force better '
        'fit; that wrecks geometry. The aim is a low-energy model that also fits.')
    return out


def set_mdff(session, model, volume=None, enabled=None, coupling_constant=None):
    '''Adjust MDFF (map-fitting) coupling state for a model.

    WARNING: raising coupling_constant to "improve fit" force-fits atoms into the
    density and typically RUINS geometry -- it is almost never the right move.
    ISOLDE's aim is to coax the model into a low-energy state that *also* fits the
    map. Use this mainly to enable/disable a map or make small, deliberate
    adjustments; prefer fixing geometry over cranking the weight.

    Args:
        model: model id/spec.
        volume: map name or id to target; if omitted, applies to all maps.
        enabled: optional bool to enable/disable coupling.
        coupling_constant: optional new global coupling weight (rarely advisable).
    '''
    s = _structure_from_spec(session, model)
    from chimerax.clipper.symmetry import get_map_mgr
    from chimerax.isolde.session_extensions import get_mdff_mgr
    mmgr = get_map_mgr(s, create=False)
    if mmgr is None:
        return {'error': 'no maps associated with this model'}
    from chimerax.isolde.cmd import agent_philosophy as ap
    results = []
    flags = []
    for v in mmgr.all_maps:
        if volume is not None and v.name != volume and ('#' + v.id_string) != str(volume):
            continue
        mgr = get_mdff_mgr(s, v, create=False)
        if mgr is None:
            continue
        prior_k = float(mgr.global_k)
        if enabled is not None:
            mgr.enabled = bool(enabled)
        if coupling_constant is not None:
            mgr.global_k = float(coupling_constant)
            # Principle 4: flag a coupling *increase* (force-fit); never block it.
            flags += ap.flag_force_fit(prior_k, float(coupling_constant))
        results.append({'volume': v.name, 'enabled': bool(mgr.enabled),
                        'coupling_constant': float(mgr.global_k)})
    out = {'model': s.atomspec, 'mdff': results}
    if flags:
        out['philosophy_flags'] = flags
    if results:
        parts = []
        if enabled is not None:
            parts.append('enabled=%s' % bool(enabled))
        if coupling_constant is not None:
            parts.append('k=%.3g' % float(coupling_constant))
        from chimerax.isolde.cmd.invoke import announce_to_user
        announce_to_user(session, '[agent] mdff %s %s (%d map%s)' % (
            s.atomspec, ' '.join(parts) or 'queried', len(results),
            '' if len(results) == 1 else 's'))
        ap.announce_high_severity(session, flags)
    return out


def bfactor_outliers(session, model, z_threshold=2.5, max_results=50):
    '''Atoms with extreme B-factors (a cheap problem sniffer, esp. post-refinement).
    Caveat: genuinely mobile sidechain atoms also show high B.'''
    import numpy
    s = _structure_from_spec(session, model)
    atoms = s.atoms
    b = atoms.bfactors
    mean = float(b.mean())
    std = float(b.std()) or 1.0
    z = (b - mean) / std
    idx = numpy.where(numpy.abs(z) >= z_threshold)[0]
    idx = idx[numpy.argsort(-numpy.abs(z[idx]))][:max_results]
    items = [{'spec': atoms[int(i)].string(style='command'),
              'bfactor': float(b[int(i)]), 'z': float(z[int(i)])} for i in idx]
    return {'model': s.atomspec, 'mean_bfactor': mean, 'std_bfactor': std,
            'z_threshold': z_threshold, 'n_outliers': int(len(idx)), 'items': items,
            'note': 'Extreme B often flags problems, but can be genuinely mobile sidechains.'}


def _atoms_from_spec(session, spec):
    from chimerax.core.commands import AtomSpecArg
    parsed, _t, _r = AtomSpecArg.parse(str(spec), session)
    return parsed.evaluate(session).atoms


def tug(session, spec, target=None, to_spec=None, spring_constant=None, release=False):
    '''
    Tug atoms toward a target using a transient spring, in a running simulation.

    Tugging is ISOLDE's gentle, geometry-preserving way to move things: it applies
    a spring pulling the selection toward a target point and lets the molecular
    dynamics relax everything around it -- the agent expresses *intent* ("this
    belongs roughly there") and MD finds a low-energy path, rather than hard-
    setting coordinates. Great for nudging a residue/fragment into density or out
    of a clash, then watching it settle. Effect plays out over sim frames.

    Requires a running simulation; the tugged atoms must be mobile (heavy atoms in
    the active sim). For a multi-atom selection the whole group is translated
    rigidly toward the target (internal geometry preserved). The spring stays
    until re-tugged or released (release=true). A very strong spring_constant
    force-fits and degrades geometry -- prefer the default and let MD do the work.

    Args:
        spec: atom-spec of the atom(s) to tug.
        target: optional [x, y, z] in Angstroms to pull the selection's centroid to.
        to_spec: alternatively, tug toward the centroid of this atom-spec.
        spring_constant: optional override (kJ/mol/A^2); default is the mouse-tug value.
        release: if true, stop tugging the selection instead of applying a tug.
    '''
    import numpy
    isolde = getattr(session, 'isolde', None)
    if isolde is None or not getattr(isolde, 'simulation_running', False):
        return {'error': 'tugging requires a running simulation; run "isolde sim start" first.'}
    sim_mgr = getattr(isolde, 'sim_manager', None)
    tugm = getattr(sim_mgr, 'tuggable_atoms_mgr', None)
    if tugm is None:
        return {'error': 'no tuggable-atoms manager (is a simulation actually running?)'}
    try:
        atoms = _atoms_from_spec(session, spec)
    except Exception as e:
        return {'error': 'bad spec %r: %s' % (spec, e)}
    if not len(atoms):
        return {'error': 'spec %r selected no atoms' % spec}

    # get_tuggables returns the collection of *tuggable* atoms in the selection
    # (heavy / sim-mobile; hydrogens and fixed atoms are dropped).
    tugs = tugm.get_tuggables(atoms)
    if not len(tugs):
        return {'error': 'none of the selected atoms are tuggable (mobile heavy '
                         'atoms in the running sim). Tug a heavy atom in the mobile selection.'}

    from chimerax.isolde.cmd.invoke import announce_to_user
    if release:
        tugs.enableds = False
        announce_to_user(session, '[agent] tug %s → released (%d)' % (spec, len(tugs)))
        return {'spec': spec, 'released': len(tugs)}

    if to_spec is not None:
        try:
            target = _atoms_from_spec(session, to_spec).coords.mean(axis=0)
        except Exception as e:
            return {'error': 'bad to_spec %r: %s' % (to_spec, e)}
    if target is None:
        return {'error': 'provide a target [x,y,z] or a to_spec to tug toward'}
    target = numpy.asarray(target, dtype=float)
    if target.shape != (3,):
        return {'error': 'target must be [x, y, z]'}

    tug_atoms = tugs.atoms
    base = tug_atoms.coords
    offset = target - base.mean(axis=0)      # rigid translation of the group
    tugs.targets = base + offset
    # Mass-weight the spring (as the mouse tug does) so all atoms respond evenly.
    # mouse_tug_spring_constant is an OpenMM Quantity (value + units); strip it to a
    # plain float so the mass-weighted assignment AND the strong-tug guardrail
    # comparison both work (float(Quantity) raises, which silently disabled the
    # guardrail).
    from openmm import unit as _unit
    default_spring = isolde.sim_params.mouse_tug_spring_constant
    if _unit.is_quantity(default_spring):
        default_spring = default_spring.value_in_unit(default_spring.unit)
    default_spring = float(default_spring)
    if spring_constant is None:
        spring_constant = default_spring
    spring_constant = float(spring_constant)
    tugs.spring_constants = spring_constant * tug_atoms.elements.masses.astype(float)
    tugs.enableds = True

    out = {'spec': spec, 'n_tugged': len(tugs), 'n_skipped': len(atoms) - len(tugs),
           'target': target.tolist(),
           'displacement_A': float(numpy.linalg.norm(offset)),
           'note': ('transient spring(s) applied; the move plays out over the next '
                    'sim frames — poll "isolde sim status" (coord_checksum) and '
                    'map_info to watch it settle. Release with tug(spec, release=true). '
                    'Tugging coaxes the model; it does not hard-set coordinates.')}
    # Principle 4/5: flag a stronger-than-default spring (force-fit risk); no block.
    from chimerax.isolde.cmd import agent_philosophy as ap
    flags = ap.flag_strong_tug(spring_constant, default_spring)
    if flags:
        out['philosophy_flags'] = flags
    announce_to_user(session, '[agent] tug %s → %.2f Å (%d atom%s)' % (
        spec, out['displacement_A'], len(tugs), '' if len(tugs) == 1 else 's'))
    ap.announce_high_severity(session, flags)
    return out


def view_state(session):
    '''Current selection + view-orientation summary (orient renders / know what is
    selected). Also reports ISOLDE's residue-stepper cursor if available.'''
    out = {}
    try:
        from chimerax.atomic import selected_atoms
        sel = selected_atoms(session)
        out['selection'] = {
            'n_atoms': len(sel),
            'residue_specs': [r.string(style='command') for r in sel.unique_residues[:25]],
        }
    except Exception:
        pass
    try:
        v = session.main_view
        out['center_of_rotation'] = [float(x) for x in v.center_of_rotation]
    except Exception:
        pass
    isolde = getattr(session, 'isolde', None)
    m = getattr(isolde, 'selected_model', None) if isolde else None
    if m is not None:
        out['isolde_selected_model'] = m.id_string
        try:
            from chimerax.isolde.navigate import get_stepper
            cur = get_stepper(m).current_residue
            out['stepper_residue'] = cur.string(style='command') if cur is not None else None
        except Exception:
            pass
    return out
