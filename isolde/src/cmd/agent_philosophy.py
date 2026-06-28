# @Author: Tristan Croll
# @Date:   28-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 28-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
The canonical ISOLDE *philosophy* for agent/MCP consumers -- one source of truth.

An ISOLDE command surface tells an agent *what it can do*; it says nothing about
*what ISOLDE is for or how it must be used*. Without that, an eyeless agent does
exactly the things ISOLDE was built to prevent: it chases validation metrics,
cranks map weights to "improve fit", treats sub-Angstrom errors as trivial, and
declares a region "done" because the numbers look good.

This module is the single home for the design philosophy (distilled from the
author's talks). It is surfaced through every layer of the agent surface, all of
which read from *here* so nothing drifts:

  - the MCP server's ``initialize`` *instructions* (always in the agent's context)
    -> :data:`INSTRUCTIONS`
  - an MCP resource ``isolde://philosophy`` (the full text, fetched on demand)
    -> :func:`principles_document`
  - the REST tool manifest preamble -> :data:`MANIFEST_PREAMBLE`
  - per-action *soft guardrails* that flag (never block) philosophy-violating moves
    -> :func:`evaluate_guardrails` and the ``flag_*`` detectors

The standalone MCP process cannot import this module (it depends only on ``mcp`` +
stdlib), so it fetches this content over REST at start-up via the
``agent_philosophy`` method (see remote_control/rest_server/agent_methods.py).
'''

# ---------------------------------------------------------------------------
# The principles (the full, canonical text -- served as the isolde://philosophy
# resource and the source of the condensed INSTRUCTIONS titles).
# ---------------------------------------------------------------------------

PRINCIPLES = [
    {'id': 1,
     'title': 'The map is the evidence; the model is a hypothesis.',
     'text': 'Your job is to make the model explain the density with good physics '
             '-- not to make the numbers look good. Treat the experimental map as '
             'ground truth and the coordinates as a falsifiable interpretation of it.'},
    {'id': 2,
     'title': '"Non-outlier" is not the same as "correct".',
     'text': 'A low clashscore, few Ramachandran/rotamer outliers, or good '
             'Rwork/Rfree do NOT prove a region is right. Validation metrics have '
             'false negatives and can be overfit or gamed. (Empirically: PDB 7adk '
             'and 8g99 are both 2.8 A with near-identical standard-metric '
             'percentiles, yet 7adk was badly wrong -- CaBLAM 14.7% disfavoured, '
             'low Q-scores, CA shifts > 2.5 A on rebuild.) A good score is never a '
             'reason to stop looking.'},
    {'id': 3,
     'title': 'Small errors are the dangerous ones.',
     'text': 'Sub-Angstrom switches, cis-peptides, flipped ligands and single bad '
             'rotamers pass static validation but have outsized, non-linear effects '
             'on dynamics, mechanism and ligand binding. ~99% of wwPDB users use '
             'structures rather than build them -- they will trust whatever is '
             'deposited. Treat small geometric details as high-stakes.'},
    {'id': 4,
     'title': 'Never force-fit.',
     'text': 'Do not raise MDFF map-coupling (or tug hard) to chase fit -- it '
             'wrecks geometry. The aim is a low-energy model that ALSO fits. If the '
             'model does not fit, fix the model, not the weight.'},
    {'id': 5,
     'title': 'Trust the physics; let it settle.',
     'text': 'ISOLDE behaves like real atoms (full electrostatics + van der Waals), '
             'which reveals otherwise-obscure detail. Express intent -- tug toward '
             'density -- and let molecular dynamics relax everything around it; do '
             'not hard-set coordinates. Settle, then check, then rebuild only where '
             'needed.'},
    {'id': 6,
     'title': 'Reference restraints supply what the density lacks -- but the reference can be wrong.',
     'text': 'Adaptive distance/torsion restraints are fuzzy "top-out": they pull '
             'only while geometry is near the reference and yield to genuine density '
             'signal. Where the experiment speaks, it wins. An AlphaFold or homologue '
             'reference may be in the wrong conformation -- do not let it override '
             'the map.'},
    {'id': 7,
     'title': 'Human eyes are the gold standard -- you work alongside, not instead.',
     'text': '"Human eyes should see each residue in density at least once." You '
             'share a live GUI with a human. Surface what you did and where you are '
             'uncertain, and never silently mark a region resolved -- the human is '
             'the safety net, so keep them able to see and veto your moves.'},
    {'id': 8,
     'title': "ISOLDE refines, but isn't the last word on fine geometry (yet).",
     'text': 'Building, real-space settling, validation, and B-factor/occupancy '
             'refinement (isolde brefine / brsr) are all in scope today. What is not '
             'there yet is reconciliation to Engh & Huber-grade bond/angle targets '
             '-- it is coming, but for now a final geometry polish in another package '
             'may still be warranted.'},
    {'id': 9,
     'title': 'Know the coverage boundaries.',
     'text': 'Real-time validation covers protein torsion-based metrics; '
             'nucleic-acid and ligand checks are partial. The force field covers '
             'protein/nucleic-acid/glycans plus ~15k common ligands and drug-like '
             'elements (C,H,N,O,S,P,F,Cl,Br,I); novel covalent linkages between '
             'residues are unsupported. Outside coverage, flag the gap -- do not '
             'fabricate confidence.'},
    {'id': 10,
     'title': 'Prefer "where do I look first?" over "is it done?"',
     'text': 'Holistic, low-false-negative triage is the goal: rank the most '
             'pressing issues (problem_zones) and focus effort there. Use the '
             'perception helpers to decide where to look, not to rubber-stamp a '
             'model as finished.'},
]


# ---------------------------------------------------------------------------
# Operational layer (distilled from ISOLDE's tutorials): the order of operations
# and rules of thumb the principles imply but don't spell out. This is "case law"
# to the principles' "constitution" -- it lives in the on-demand resource/document,
# NOT in the always-in-context INSTRUCTIONS.
# ---------------------------------------------------------------------------

WORKFLOW = [
    'Prepare for simulation: add hydrogens (addh), remove alternate conformations, '
    'and resolve unparametrised residues (isolde preflight) so every residue is '
    'known to the force field.',
    'Settle once. Run a whole-model simulation to relieve clashes and let obvious '
    'errors fix themselves. You typically simulate the whole model only twice -- '
    'once at the start to relieve clashes, once at the end to settle.',
    'Triage with problem_zones: find the worst clusters and fix the biggest first -- '
    'many downstream issues cascade away once the root cause is fixed.',
    'Fix locally: peptide flips, rotamer corrections, register shifts. Tug in short '
    'bursts to *help* atoms into density, then let MD settle -- never hard-pull.',
    'Reconcile restraints with the map. Where reference/AlphaFold restraints fight '
    'strong density (a cluster of strained restraints plus poor fit), the reference '
    'is probably wrong: release those restraints and let the density win.',
    'Final settle: a whole-model simulation with temperature dropped toward zero, to '
    'bring local geometry to equilibrium before writing coordinates.',
    'Inspect everything: step through every residue in context with its density '
    '(isolde stepto / navigation). Your eyes remain the gold standard.',
    'Refine externally: brefine/brsr handle B-factors and occupancies; a final '
    'fine-geometry polish in a dedicated package may still be warranted.',
]

RULES_OF_THUMB = [
    'Never add or remove atoms while a simulation is running.',
    'Do not chase fit by raising MDFF coupling. Watch for over-fitting symptoms: '
    'backbone twisting out of stable secondary structure, or sidechains pulled too '
    'aggressively into the map.',
    'Tug in short bursts to help the model into place; one hard pull force-fits and '
    'wrecks local geometry.',
    'Severe clashes (atoms passing through each other) will not resolve by naive '
    'minimisation. Separate them gently -- soft-core nonbonded potentials plus local '
    'restraints over a region big enough for the pieces to slide apart -- rather than '
    'cranking forces.',
    'Validation is a hint, not a verdict: not every outlier is wrong, and not being '
    'an outlier does not make something correct. A peptide twist beyond ~30 deg is '
    'effectively never real; a non-proline cis bond is rare (~3 in 10,000) and, when '
    'genuine, well-resolved and functionally interesting; ~5% of prolines are '
    'legitimately cis. Always confirm against the density.',
    'AlphaFold / reference models: residues with pLDDT < 50 are essentially junk -- do '
    'not interpret their coordinates; only trust or restrain residue pairs with low '
    'predicted error (PAE no more than ~4 Angstrom). A confident prediction can still '
    'be wrong where the map disagrees.',
    'Out-of-register stretches: use the register-shifter, not a tug-of-war between '
    'distance restraints.',
    'For final geometry, settle at low or zero temperature.',
]


# ---------------------------------------------------------------------------
# The condensed "constitution" -- always in the agent's context (MCP instructions).
# Built from the principle TITLES so it can never drift from PRINCIPLES; only the
# framing prose is hand-written. The full text is the isolde://philosophy resource.
# ---------------------------------------------------------------------------

_INSTRUCTIONS_PREAMBLE = (
    'You are driving ISOLDE: interactive, physics-based (all-atom MD) model '
    'building and refinement into crystallographic or cryo-EM density, alongside '
    'a human in a live GUI. Before you act, internalise these principles -- they '
    'describe the nature of the target and how to avoid getting the user into a '
    'mess. They matter more than any single tool call.')

_INSTRUCTIONS_CODA = (
    'Mutating actions and their outcomes are echoed to the human\'s GUI log; '
    'read-only queries are not. Some results carry "philosophy_flags" -- advisory '
    'notes (never blocks) when a move looks like it violates a principle above; '
    'heed them. The isolde://philosophy resource holds the full text plus a '
    'recommended end-to-end workflow and rules of thumb -- read it before planning '
    'a build.')


def _build_instructions():
    lines = [_INSTRUCTIONS_PREAMBLE, '']
    for p in PRINCIPLES:
        lines.append('%2d. %s' % (p['id'], p['title']))
    lines += ['', _INSTRUCTIONS_CODA]
    return '\n'.join(lines)


INSTRUCTIONS = _build_instructions()

MANIFEST_PREAMBLE = (
    'These tools drive ISOLDE (physics-based model building/refinement into '
    'density, with a human in a live GUI). Before mutating, read the principles: '
    'the map is the evidence, not the scores; never force-fit; small errors '
    'matter; let the physics settle; you work alongside a human. Full text: the '
    'isolde://philosophy MCP resource (or the agent_philosophy REST method).')


def principles_document():
    '''The full philosophy as a single text block (the isolde://philosophy
    resource body): the principles, then the recommended workflow and rules of
    thumb (the operational "case law").'''
    lines = ['ISOLDE -- design philosophy for agents', '']
    for p in PRINCIPLES:
        lines.append('%d. %s' % (p['id'], p['title']))
        lines.append('   %s' % p['text'])
        lines.append('')
    lines.append('Recommended workflow')
    lines.append('')
    for i, step in enumerate(WORKFLOW, 1):
        lines.append('%d. %s' % (i, step))
    lines.append('')
    lines.append('Rules of thumb')
    lines.append('')
    for rule in RULES_OF_THUMB:
        lines.append('- %s' % rule)
    return '\n'.join(lines).rstrip() + '\n'


def as_payload():
    '''The dict returned by the agent_philosophy REST method (the MCP server's
    single fetch point).'''
    return {
        'instructions': INSTRUCTIONS,
        'preamble': MANIFEST_PREAMBLE,
        'principles': [dict(p) for p in PRINCIPLES],
        'workflow': list(WORKFLOW),
        'rules_of_thumb': list(RULES_OF_THUMB),
        'document': principles_document(),
    }


# ---------------------------------------------------------------------------
# Soft guardrails: detect a philosophy-violating move and return advisory flags.
# They NEVER block. Messages are constants here (tagged with a principle id) so
# the scattered per-tool warnings can converge on one source and not drift.
# ---------------------------------------------------------------------------

# Tunable threshold: a single command shifting a heavy atom more than this (while
# not just letting a sim settle) is worth a "did you mean to do that? confirm
# against density / let the human see it" nudge. First guess; verify live.
LARGE_MOVE_ANGSTROM = 5.0

_MSG = {
    'force_fit': ('Raised MDFF map-coupling. Principle 4: raising coupling to chase '
                  'fit force-fits atoms and wrecks geometry -- fix the model, not the '
                  'weight. The aim is a low-energy model that also fits.'),
    'strong_tug': ('Tug spring stronger than the default. Principle 4/5: a strong tug '
                   'force-fits and degrades local geometry -- prefer the default and '
                   'let MD relax the surroundings.'),
    'large_move': ('Large single-step atom shift (%.1f A). Principle 5/7: prefer '
                   'gentle settling over big shoves, confirm the move against the '
                   'density, and make sure the human can see it.'),
    'clean_metric': ('Validation reports no outliers here. Principle 2: "non-outlier" '
                     'is not "correct" -- metrics have false negatives and can be '
                     'overfit. A clean score is not a reason to stop looking; confirm '
                     'against the map.'),
}


def _flag(principle, severity, message):
    return {'principle': principle, 'severity': severity, 'message': message}


def flag_force_fit(prior_k, new_k):
    '''High-severity flag if an MDFF coupling change raises the weight.'''
    try:
        if prior_k is not None and new_k is not None and float(new_k) > float(prior_k):
            return [_flag(4, 'high', _MSG['force_fit'])]
    except (TypeError, ValueError):
        pass
    return []


def flag_strong_tug(spring_constant, default):
    '''Warn if a tug uses a spring stronger than the default mouse-tug value.'''
    try:
        if spring_constant is not None and default is not None \
                and float(spring_constant) > float(default):
            return [_flag(4, 'warn', _MSG['strong_tug'])]
    except (TypeError, ValueError):
        pass
    return []


def flag_large_move(witness):
    '''Warn if a coord_move witness recorded a large single-step heavy-atom shift.'''
    if not isinstance(witness, dict):
        return []
    shift = witness.get('max_atom_shift')
    try:
        if shift is not None and float(shift) > LARGE_MOVE_ANGSTROM:
            return [_flag(5, 'warn', _MSG['large_move'] % float(shift))]
    except (TypeError, ValueError):
        pass
    return []


# Validation result keys that count genuine problems (tolerant across report shapes).
_OUTLIER_COUNT_KEYS = ('n_outlier', 'n_outliers', 'n_clashes', 'n_clash', 'n_twisted')


def flag_clean_validation(result):
    '''Info-level caveat when a validation query reports zero problems: a clean
    score is never proof of correctness (principle 2).'''
    if not isinstance(result, dict):
        return []
    counts = [result[k] for k in _OUTLIER_COUNT_KEYS
              if isinstance(result.get(k), (int, float))]
    if counts and sum(counts) == 0:
        return [_flag(2, 'info', _MSG['clean_metric'])]
    return []


def evaluate_guardrails(category=None, witness=None, result=None):
    '''Guardrails evaluable at the universal invoke choke point, from the command
    category + its witness + its (jsonified) result. The force-fit and strong-tug
    detectors need the specific prior/args and are called directly by their helpers
    (set_mdff/tug); they are not reachable from here.'''
    flags = []
    flags += flag_large_move(witness)
    if category == 'validation':
        flags += flag_clean_validation(result)
    return flags


def announce_high_severity(session, flags):
    '''Echo any high-severity flags to the human's GUI log (principle 7: keep the
    human able to see and veto). Best-effort; never raises into the command path.'''
    if not flags:
        return
    try:
        from .invoke import announce_to_user
    except Exception:
        return
    for f in flags:
        if f.get('severity') == 'high':
            announce_to_user(session, '[agent] ⚠ %s' % f.get('message', ''))
