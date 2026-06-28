# @Author: Tristan Croll
# @Date:   27-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 27-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Opt-in allowlist of ISOLDE commands exposed to agents (design decision 3).

Auto-discovery (see :mod:`.command_registry`) keeps the command *manifest*
current with zero effort, but **exposure to agents is opt-in**: a command is
agent-reachable only when it is explicitly listed here with ``agent_safe=True``.
A newly added or half-baked command -- including ones that arrive on a future
branch -- is captured by the registry but stays invisible to agents until a
human adds a line below. That is the smallest footgun surface and matches the
security posture around ChimeraX's arbitrary-code command interface.

Each entry maps a full command name to a dict with optional keys:
    category      -- grouping shown in the manifest, e.g. 'simulation'
    agent_safe    -- bool; only True commands are reachable (default False)
    witness       -- post-condition witness id (see .witness), e.g. 'coord_move'
    long_running  -- bool; routed through the async job pattern (jobs.py)

Entries for commands that are not (yet) registered are harmless: they sit
dormant until the matching command appears, so manifests for unmerged branches
(e.g. rdkit's chirality / ligand commands) can be pre-staged here.
'''

# Category-level defaults. A command inherits its category's flags unless it
# sets its own. Keeping witness/long_running at the category level avoids
# bolting metadata onto individual commands where the whole class behaves alike.
CATEGORY_DEFAULTS = {
    'simulation':   dict(agent_safe=True),
    'query':        dict(agent_safe=True),          # read-only; no witness needed
    'validation':   dict(agent_safe=True),          # read-only structured reports
    'manipulation': dict(agent_safe=True, witness='coord_move'),
    'building':     dict(agent_safe=True, witness='coord_move'),
    # Restraint commands change restraint STATE (not coordinates directly) and
    # mostly return None -> the restraint_change witness counts restraints added/
    # released/enabled via the managers (the silent-no-op guard for restraints).
    'restraints':   dict(agent_safe=True, witness='restraint_change'),
    'navigation':   dict(agent_safe=True),
    'refinement':   dict(agent_safe=True, witness='convergence', long_running=True),
    'parameterise': dict(agent_safe=True),
    'output':       dict(agent_safe=True),
}


# Per-command opt-in. `category` selects the defaults above; any of agent_safe /
# witness / long_running given here override the category default for that one
# command. NOTE: this is the ONLY place exposure is granted -- adding a command
# to the registry does not expose it.
AGENT_COMMANDS = {
    # --- Simulation control ---
    'isolde start':              dict(category='simulation'),
    'isolde sim':                dict(category='simulation', witness='sim_control'),
    'isolde set':                dict(category='simulation'),
    'isolde select':             dict(category='simulation'),
    'isolde report':             dict(category='query'),
    'isolde ignore':             dict(category='simulation'),
    'isolde ~ignore':            dict(category='simulation'),

    # --- Query / validation (query plane) ---
    'isolde status':             dict(category='query'),
    'isolde sim status':         dict(category='query'),
    'isolde validate peptidebonds': dict(category='validation'),
    'isolde validate rama':      dict(category='validation'),
    'isolde validate rotamers':  dict(category='validation'),
    'isolde validate clashes':   dict(category='validation'),
    'isolde validate chirals':   dict(category='validation'),
    # Read-only MD-readiness checks (run before starting a simulation).
    'isolde preflight hydrogens':  dict(category='query'),
    'isolde preflight parameters': dict(category='query'),
    'isolde preflight disulfides': dict(category='query'),
    'isolde preflight altlocs':    dict(category='query'),

    # --- Restraints (guided model building) ---
    'isolde restrain distances':       dict(category='restraints'),
    'isolde restrain single distance': dict(category='restraints'),
    'isolde restrain torsions':        dict(category='restraints'),
    'isolde restrain ss':              dict(category='restraints'),
    'isolde restrain basepairs':       dict(category='restraints'),
    'isolde restrain ligands':         dict(category='restraints'),
    'isolde adjust distances':         dict(category='restraints'),
    'isolde adjust torsions':          dict(category='restraints'),
    'isolde release distances':        dict(category='restraints'),
    'isolde release torsions':         dict(category='restraints'),

    # --- B-factor / occupancy refinement (long-running; async job pattern) ---
    'isolde brefine':                dict(category='refinement'),
    'isolde brefine optimiseparams': dict(category='refinement'),
    'isolde brsr':                   dict(category='refinement'),

    # --- Structure building / modification ---
    'isolde pepflip':            dict(category='manipulation'),
    'isolde cisflip':            dict(category='manipulation'),
    'isolde modify his':         dict(category='manipulation'),
    'isolde add water':          dict(category='building'),
    'isolde add aa':             dict(category='building'),
    'isolde adjust bfactors':    dict(category='building'),
    'isolde add disulfides auto': dict(category='building'),
    'isolde clear altlocs':      dict(category='building'),

    # --- Navigation / view ---
    'isolde stepto':             dict(category='navigation'),
    'isolde jumpto':             dict(category='navigation'),

    # --- Chirality (rdkit branch) ---
    # Correction tool; the matching read-only validator is
    # 'isolde validate chirals' (in the query/validation block above). The live
    # markup command ('isolde annotate chirals', like the other 'isolde annotate'
    # markup) is GUI-only and intentionally not exposed to agents.
    'isolde chiralflip':         dict(category='manipulation'),
}


def metadata_for(name: str) -> dict:
    '''
    Resolve the effective agent metadata for command *name*: category defaults
    merged with any per-command overrides. Returns an empty-ish dict (agent_safe
    False) for commands not listed here -- i.e. captured but not exposed.
    '''
    entry = AGENT_COMMANDS.get(name)
    if entry is None:
        return {'agent_safe': False}
    category = entry.get('category')
    merged = dict(CATEGORY_DEFAULTS.get(category, {}))
    merged.update({k: v for k, v in entry.items() if k != 'category'})
    merged['category'] = category
    merged.setdefault('agent_safe', False)
    return merged
