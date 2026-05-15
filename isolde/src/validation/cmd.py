# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 14-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll



def rota(session, structures=None, report=False):
    '''
    Add a live rotamer validator to each of the given structures, and optionally
    report a summary of current outliers.
    '''
    from chimerax.isolde import session_extensions as sx
    if structures is None:
        from chimerax.atomic import AtomicStructure
        structures = [m for m in session.models.list() if type(m)==AtomicStructure]
    for structure in structures:
        sx.get_rota_annotator(structure)
    if report:
        report_str = 'NON-FAVOURED ROTAMERS: \n'
        for structure in structures:
            data = _compute_rotamer_report(session, structure,
                include='nonfavored')
            for it in data['items']:
                report_str += '#{:<6} {}:\t{} {} (P={:.4f})\n'.format(
                    structure.id_string, it['chain_id'], it['name'],
                    it['number'], it['score']
                )
        session.logger.info(report_str)

def unrota(session, structures=None):
    '''
    Delete any rotamer annotators associated with the given models.
    '''
    if structures is None:
        from chimerax.atomic import AtomicStructure
        structures = [m for m in session.models.list() if type(m)==AtomicStructure]
    from chimerax.isolde import session_extensions as sx
    for structure in structures:
        ra = sx.get_rota_annotator(structure, create=False)
        if ra is not None:
            session.models.close([ra])

def rama(session, structures=None, show_favored=True, report=False):
    '''
    Add a live Ramachandran validator to each of the given structures, and
    optionally report a summary of current outliers and cis/twisted peptide
    bonds.
    '''
    from chimerax.isolde import session_extensions as sx
    if structures is None:
        from chimerax.atomic import AtomicStructure
        structures = [m for m in session.models.list() if type(m)==AtomicStructure]
    for structure in structures:
        ra = sx.get_rama_annotator(structure)
        ra.hide_favored = not show_favored
    if report:
        pep_data_by_structure = [
            (structure, _compute_peptide_bond_report(session, structure))
            for structure in structures
        ]
        report_str = 'RAMACHANDRAN OUTLIERS: \n'
        for structure in structures:
            rd = _compute_rama_report(session, structure, include='outliers')
            for it in rd['items']:
                report_str +='#{:<6} {}:\t{} {}\n'.format(
                    structure.id_string, it['chain_id'],
                    it['name'], it['number'])
        report_str += '\nCIS PEPTIDE BONDS: \n'
        for structure, pd in pep_data_by_structure:
            for it in pd['items']:
                if it['conformation'] != 'cis':
                    continue
                r2 = it['res2']
                report_str +='#{:<6} {}:\t{} {}\n'.format(
                    structure.id_string, r2['chain_id'],
                    r2['name'], r2['number']
                )
        report_str += '\nTWISTED PEPTIDE BONDS: \n'
        for structure, pd in pep_data_by_structure:
            for it in pd['items']:
                if it['conformation'] != 'twisted':
                    continue
                r2 = it['res2']
                report_str += '#{:<6} {}:\t{} {} ({:.1f}°)\n'.format(
                    structure.id_string, r2['chain_id'],
                    r2['name'], r2['number'], it['omega_deg']
                )
        session.logger.info(report_str)

def unrama(session, structures=None):
    '''
    Delete any Ramachandran annotators associated with the given models.
    '''
    if structures is None:
        from chimerax.atomic import AtomicStructure
        structures = [m for m in session.models.list() if type(m)==AtomicStructure]
    from chimerax.isolde import session_extensions as sx
    for structure in structures:
        ra = sx.get_rama_annotator(structure, create=False)
        if ra is not None:
            session.models.close([ra])

def register_rota(logger):
    from chimerax.core.commands import (
        register, CmdDesc, BoolArg, create_alias
    )
    from chimerax.atomic import StructuresArg
    desc = CmdDesc(
        optional=[
            ('structures', StructuresArg),
            ],
        keyword=[
            ('report', BoolArg),
        ],
        synopsis='Add rotamer validator markup to models and optionally report current outliers'
    )
    register('rota', desc, rota, logger=logger)
    undesc = CmdDesc(
        optional=[('structures', StructuresArg)],
        synopsis='Close the rotamer annotators for the given models (or all if no models given)'
    )
    register('rota stop', undesc, unrota, logger=logger)
    create_alias('~rota', 'rota stop $*', logger=logger)


def register_rama(logger):
    from chimerax.core.commands import (
        register, CmdDesc, BoolArg, create_alias
    )
    from chimerax.atomic import StructuresArg
    desc = CmdDesc(
        optional=[
            ('structures', StructuresArg),
            ],
        keyword=[
            ('show_favored', BoolArg),
            ('report', BoolArg),
        ],
        synopsis='Add Ramachandran validator markup to models and optionally report current outliers'
    )
    register('rama', desc, rama, logger=logger)
    undesc = CmdDesc(
        optional=[('structures', StructuresArg)],
        synopsis='Close the Ramachandran annotators for the given models (or all if no models given)'
    )
    register('rama stop', undesc, unrama, logger=logger)
    create_alias('~rama', 'rama stop $*', logger=logger)


# ---------------------------------------------------------------------------
# Read-only ``isolde preflight`` commands.
#
# These commands let an agent (or scripted user) ask whether a model is ready
# for an ISOLDE simulation without actually trying to start one. They are safe
# to call at any time: they never construct an OpenMM ``Context``, never raise
# on missing parameters, and never modify the model.
# ---------------------------------------------------------------------------

from chimerax.core.errors import UserError


# Subcommand names exposed under ``isolde preflight``. Used both for
# registration order and to build the helpful error raised by the bare
# ``isolde preflight`` parent handler below.
_PREFLIGHT_SUBCOMMANDS = ('hydrogens', 'parameters', 'disulfides', 'altlocs')


def isolde_preflight(session, model=None):
    '''
    Bare ``isolde preflight`` handler. Always raises ``UserError`` listing
    the available subcommands.

    Registered as the parent of ``isolde preflight hydrogens`` /
    ``parameters`` / ``disulfides`` / ``altlocs`` so that calls like
    ``isolde preflight #1.2`` (model spec but no subcommand) get a useful
    "expected one of these" message instead of ChimeraX's generic
    ``Unknown command: isolde preflight #1.2``.
    '''
    raise UserError(
        "'isolde preflight' requires a subcommand. Available: "
        + ', '.join(_PREFLIGHT_SUBCOMMANDS)
        + ". Example: 'isolde preflight hydrogens #1'."
    )

from .unparameterised import (
    H_TO_HEAVY_ATOM_THRESHOLD_RATIO,
    suspiciously_low_h,
    waters_without_h,
)


def _resolve_model(session, model):
    if model is not None:
        return model
    isolde = getattr(session, 'isolde', None)
    if isolde is not None and isolde.selected_model is not None:
        return isolde.selected_model
    raise UserError(
        'No model specified and no model is currently selected in ISOLDE. '
        'Either pass a model argument or run "isolde select <model>" first.'
    )


def _residue_summary(r):
    return {
        'chain_id': r.chain_id,
        'name': r.name,
        'number': int(r.number),
        'insertion_code': r.insertion_code,
        'spec': r.atomspec,
    }


def _residue_label(r):
    icode = r.insertion_code.strip()
    return '/{} {} {}{}'.format(r.chain_id, r.name, r.number, icode)


def isolde_preflight_hydrogens(session, model=None):
    '''
    Preflight check: does ``model`` (or ISOLDE's currently selected model)
    appear to have a complete set of hydrogens for an MD simulation? Uses the
    same heuristics ISOLDE's "Unparametrised residues" panel applies before
    launching a simulation.

    Does not modify the model and does not start a simulation. Returns a dict
    suitable for programmatic consumption; also logs a one-line summary.
    '''
    m = _resolve_model(session, model)
    residues = m.residues
    atoms = residues.atoms
    enames = atoms.element_names
    n_h = int((enames == 'H').sum())
    n_heavy = int((enames != 'H').sum())
    ratio = (n_h / n_heavy) if n_heavy else 0.0
    waters = residues[residues.names == 'HOH']
    waters_missing = sum(1 for w in waters if len(w.atoms) != 3)

    if n_h == 0:
        status = 'missing'
        message = 'No hydrogens present.'
    elif suspiciously_low_h(residues):
        message = (
            'H/heavy-atom ratio {:.2f} is below the expected threshold of {:.2f}.'
            .format(ratio, H_TO_HEAVY_ATOM_THRESHOLD_RATIO)
        )
        status = 'low'
    elif waters_without_h(residues):
        status = 'waters_missing'
        message = '{} water residue(s) are missing one or both hydrogens.'.format(
            waters_missing
        )
    else:
        status = 'ok'
        message = 'Hydrogens look complete.'

    log = session.logger
    summary = (
        'ISOLDE hydrogen check ({}): {} ({} H / {} heavy, ratio {:.2f}; '
        '{} waters, {} missing H).'.format(
            m.atomspec, message, n_h, n_heavy, ratio, len(waters), waters_missing
        )
    )
    if status == 'ok':
        log.info(summary)
    else:
        log.warning(summary)
        log.info('Recommendation: run "addh" before starting an ISOLDE simulation.')

    return {
        'model': m.atomspec,
        'status': status,
        'message': message,
        'hydrogen_count': n_h,
        'heavy_atom_count': n_heavy,
        'ratio': ratio,
        'water_count': int(len(waters)),
        'waters_missing_hydrogens': int(waters_missing),
        'recommend_addh': status != 'ok',
    }


def _get_forcefield(session, ff_name):
    '''
    Return ``(forcefield, ligand_db, ff_name)``. Reuses ISOLDE's cached
    ``ForcefieldMgr`` if available so we don't reload the forcefield from
    disk; otherwise creates a temporary local manager. Either way no
    simulation is started.
    '''
    isolde = getattr(session, 'isolde', None)
    if isolde is not None:
        ffmgr = isolde.forcefield_mgr
        if ff_name is None:
            ff_name = isolde.sim_params.forcefield
    else:
        from ..openmm.forcefields import ForcefieldMgr
        ffmgr = ForcefieldMgr(session)
        if ff_name is None:
            from ..constants import defaults
            ff_name = defaults.OPENMM_FORCEFIELD
    ff = ffmgr[ff_name]
    ligand_db = ffmgr.ligand_db(ff_name)
    return ff, ligand_db, ff_name


def isolde_preflight_parameters(session, model=None, forcefield=None,
        ignore_external_bonds=True):
    '''
    Preflight check: run ISOLDE's MD-template assignment over the given model
    (or ISOLDE's currently selected model) without starting a simulation, and
    report any residues that are unparametrised or ambiguous.

    This is the same dry-run that the "Unparametrised residues" panel
    performs - it does not call the OpenMM ``Context`` constructor and so
    cannot raise the cryptic errors that arise from a real simulation start.

    Parameters
    ----------
    model : AtomicStructure, optional
        The model to check. Defaults to ISOLDE's currently selected model.
    forcefield : str, optional
        The forcefield name to check against (e.g. ``'amber14'``). Defaults
        to ISOLDE's currently configured forcefield.
    ignore_external_bonds : bool, default True
        Whether to ignore external bonds when matching residues to templates.
        ``True`` matches the behaviour used by the GUI panel.

    Returns
    -------
    dict
        Summary including counts and per-residue details for any unmatched or
        ambiguous residues.
    '''
    m = _resolve_model(session, model)
    log = session.logger

    ff, ligand_db, ff_name = _get_forcefield(session, forcefield)

    from chimerax.atomic import Residues
    residues = Residues(sorted(m.residues,
        key=lambda r: (r.chain_id, r.number, r.insertion_code)))

    from ..openmm.openmm_interface import (
        find_residue_templates, create_openmm_topology
    )
    template_dict = find_residue_templates(
        residues, ff, ligand_db=ligand_db, logger=log
    )
    top, residue_templates = create_openmm_topology(
        residues.atoms, template_dict
    )
    unique, ambiguous, unmatched = ff.assignTemplates(
        top, ignoreExternalBonds=ignore_external_bonds,
        explicit_templates=residue_templates,
    )

    unmatched_info = []
    for r in unmatched:
        cx_res = residues[r.index]
        info = _residue_summary(cx_res)
        # Suggest possible templates so the agent has something actionable.
        try:
            by_name, by_comp = ff.find_possible_templates(r)
        except Exception:
            by_name, by_comp = [], []
        info['candidates_by_name'] = [
            {'template': t, 'score': float(s)} for (t, s) in by_name
        ]
        info['candidates_by_topology'] = [
            {'template': t, 'score': float(s)} for (t, s) in by_comp
        ]
        unmatched_info.append(info)

    ambiguous_info = []
    for r, template_info in ambiguous.items():
        cx_res = residues[r.index]
        info = _residue_summary(cx_res)
        info['candidates'] = [t[0].name for t in template_info]
        ambiguous_info.append(info)

    n_total = len(residues)
    n_unique = len(unique)
    n_ambig = len(ambiguous)
    n_unmatched = len(unmatched)

    header = (
        'ISOLDE parameter check ({}, forcefield {}): '
        '{} residues, {} matched, {} ambiguous, {} unmatched.'.format(
            m.atomspec, ff_name, n_total, n_unique, n_ambig, n_unmatched
        )
    )
    if n_unmatched == 0 and n_ambig == 0:
        log.info(header)
    else:
        log.warning(header)
        if unmatched_info:
            log.info('  Unmatched: ' + ', '.join(
                _residue_label(residues[u.index]) for u in unmatched
            ))
        if ambiguous_info:
            log.info('  Ambiguous: ' + ', '.join(
                _residue_label(residues[r.index]) for r in ambiguous
            ))

    return {
        'model': m.atomspec,
        'forcefield': ff_name,
        'ignore_external_bonds': bool(ignore_external_bonds),
        'n_total': int(n_total),
        'n_matched': int(n_unique),
        'n_ambiguous': int(n_ambig),
        'n_unmatched': int(n_unmatched),
        'ready_for_simulation': n_unmatched == 0 and n_ambig == 0,
        'unmatched': unmatched_info,
        'ambiguous': ambiguous_info,
    }


def isolde_preflight_disulfides(session, model=None):
    '''
    Preflight check: does ``model`` (or ISOLDE's currently selected model)
    contain pairs of cysteines whose SG atoms are close enough to be
    disulfide-bonded but for which no SG-SG bond is currently present in the
    model?

    This is the same geometric check that fires the "create disulfides?" GUI
    popup the first time a model is selected in ISOLDE. Calling this command
    does *not* create any bonds, but it stamps a per-model flag so that the
    auto-popup will not subsequently fire for the same model — i.e. running
    this preflight is treated as the user/agent's acknowledgement of the
    situation. Pair this with ``isolde add disulfides auto`` if you want the
    bonds created.

    Returns a dict with three lists (``current``, ``possible``, ``ambiguous``).
    Each list entry is itself a list of residue summaries.
    '''
    m = _resolve_model(session, model)
    log = session.logger

    from ..atomic.building.build_utils import current_and_possible_disulfides
    current, possible, ambiguous = current_and_possible_disulfides(
        m, cutoff_distance=2.3
    )

    def _pair(rset):
        return [_residue_summary(r) for r in sorted(
            rset, key=lambda r: (r.chain_id, r.number, r.insertion_code)
        )]

    current_info = [_pair(p) for p in current]
    possible_info = [_pair(p) for p in possible]
    ambiguous_info = [_pair(p) for p in ambiguous]

    def _label_set(rset):
        return '-'.join(_residue_label(r) for r in sorted(
            rset, key=lambda r: (r.chain_id, r.number, r.insertion_code)
        ))

    summary = (
        'ISOLDE disulfide check ({}): {} existing, {} possible new, '
        '{} ambiguous cluster(s).'.format(
            m.atomspec, len(current), len(possible), len(ambiguous)
        )
    )
    if not possible and not ambiguous:
        log.info(summary)
    else:
        log.warning(summary)
        if possible:
            log.info('  Possible: ' + ', '.join(
                _label_set(p) for p in possible
            ))
            log.info(
                '  To create them, run "isolde add disulfides auto {}".'
                .format(m.atomspec)
            )
        if ambiguous:
            log.info('  Ambiguous (3+ cysteines clustered, manual fix required): '
                + '; '.join(_label_set(p) for p in ambiguous))

    # Calling the preflight counts as acknowledging the situation; this
    # suppresses the GUI popup that would otherwise fire on the next frame
    # after `isolde select`.
    m._isolde_disulfide_check_done = True

    return {
        'model': m.atomspec,
        'current': current_info,
        'possible': possible_info,
        'ambiguous': ambiguous_info,
        'n_current': len(current),
        'n_possible': len(possible),
        'n_ambiguous': len(ambiguous),
        'recommend_create': len(possible) > 0,
    }


def isolde_preflight_altlocs(session, model=None):
    '''
    Preflight check: does ``model`` (or ISOLDE's currently selected model)
    contain atoms with alternate conformations? ISOLDE cannot see alt locs
    during a simulation, but they are carried through to the output, so in
    most refinement workflows they should be removed before starting.

    This is the situation that triggers the "remove alt locs?" GUI popup
    inside ``Isolde.selected_model`` the first time a model is selected.
    Calling this preflight stamps a per-model flag so that the auto-popup
    will not fire for the model. Pair with ``isolde clear altlocs`` to
    actually drop them.

    Returns a dict with the per-residue list and counts.
    '''
    m = _resolve_model(session, model)
    log = session.logger

    atoms_with_altlocs = m.atoms[m.atoms.num_alt_locs > 0]
    n_atoms = int(len(atoms_with_altlocs))
    affected_residues = atoms_with_altlocs.unique_residues
    residues_info = [_residue_summary(r) for r in affected_residues]

    summary = (
        'ISOLDE altloc check ({}): {} atom(s) with alternate conformers '
        'across {} residue(s).'.format(m.atomspec, n_atoms, len(residues_info))
    )
    if n_atoms == 0:
        log.info(summary)
    else:
        log.warning(summary)
        log.info('  Affected residues: ' + ', '.join(
            _residue_label(r) for r in affected_residues
        ))
        log.info(
            '  To drop alt locs and reset occupancies, run '
            '"isolde clear altlocs {}".'.format(m.atomspec)
        )

    # See note in isolde_preflight_disulfides: the call itself acknowledges
    # the situation and suppresses the auto-popup.
    m._isolde_altloc_check_done = True

    return {
        'model': m.atomspec,
        'atoms_with_altlocs': n_atoms,
        'residues': residues_info,
        'n_residues': len(residues_info),
        'recommend_clear': n_atoms > 0,
    }


def register_preflight_commands(logger):
    from chimerax.core.commands import (
        register, CmdDesc, BoolArg, StringArg,
    )
    from ..cmd.argspec import IsoldeStructureArg

    desc_h = CmdDesc(
        optional=[('model', IsoldeStructureArg)],
        synopsis='Preflight check: does the model have plausible hydrogens for MD?',
    )
    register('isolde preflight hydrogens', desc_h,
        isolde_preflight_hydrogens, logger=logger)

    desc_p = CmdDesc(
        optional=[('model', IsoldeStructureArg)],
        keyword=[
            ('forcefield', StringArg),
            ('ignore_external_bonds', BoolArg),
        ],
        synopsis='Preflight check: do all residues match an MD forcefield template?',
    )
    register('isolde preflight parameters', desc_p,
        isolde_preflight_parameters, logger=logger)

    desc_d = CmdDesc(
        optional=[('model', IsoldeStructureArg)],
        synopsis=('Preflight check: are there cysteine pairs likely to need '
            'disulfide bonds? Suppresses the corresponding GUI popup.'),
    )
    register('isolde preflight disulfides', desc_d,
        isolde_preflight_disulfides, logger=logger)

    desc_a = CmdDesc(
        optional=[('model', IsoldeStructureArg)],
        synopsis=('Preflight check: does the model contain alternate '
            'conformers? Suppresses the corresponding GUI popup.'),
    )
    register('isolde preflight altlocs', desc_a,
        isolde_preflight_altlocs, logger=logger)

    # Parent command: catches ``isolde preflight`` (no subcommand) and
    # ``isolde preflight <model>`` (model spec but no subcommand) and
    # turns them into a helpful "expected one of: ..." error instead of
    # ChimeraX's generic ``Unknown command``. The ``optional`` model arg
    # exists so the parser cleanly consumes the spec - the handler always
    # raises before doing any work.
    desc_top = CmdDesc(
        optional=[('model', IsoldeStructureArg)],
        synopsis='Run an ISOLDE preflight check (requires a subcommand: {}).'.format(
            ', '.join(_PREFLIGHT_SUBCOMMANDS)),
    )
    register('isolde preflight', desc_top, isolde_preflight, logger=logger)


# ---------------------------------------------------------------------------
# Read-only ``isolde validate`` commands.
#
# These commands run the same scoring/validators as ISOLDE's GUI Validate
# tab and return structured results. They never modify the model and never
# start a simulation - safe to call at any time once a model is selected.
#
# Each handler returns a ``dict`` with summary counts plus a per-item list
# (suitable for direct consumption by an agent / MCP caller). Output is
# also routed through three shared keywords:
#
# - ``log`` (bool, default False): dump the full per-item table to the
#   ChimeraX HTML log wrapped in ``<pre>...</pre>`` (same pattern as the
#   ChimeraX ``clashes`` / ``hbonds`` commands).
# - ``save_file`` (path, default None): write the full table to disk.
#   ``.json`` paths get a structured JSON dump; anything else gets a
#   plain-text aligned table written via ``chimerax.io.open_output``.
# - ``limit`` (int, default None): cap the per-item list in the returned
#   dict so a giant structure doesn't blow up the agent's context. The
#   file output ignores this limit and always contains the full list.
# ---------------------------------------------------------------------------


# Don't dump arbitrarily large tables into the HTML log; clip to this many
# rows when ``log=True`` so the Reply Log stays usable.
_LOG_TABLE_ROW_LIMIT = 500


# Subcommand names exposed under ``isolde validate``. Used both for
# registration order and to build the helpful error raised by the bare
# ``isolde validate`` parent handler below.
_VALIDATE_SUBCOMMANDS = ('peptidebonds', 'rama', 'rotamers', 'clashes')


def isolde_validate(session, model=None):
    '''
    Bare ``isolde validate`` handler. Always raises ``UserError`` listing
    the available subcommands.

    Registered as the parent of ``isolde validate peptidebonds`` /
    ``rama`` / ``rotamers`` / ``clashes`` so that calls like
    ``isolde validate #1.2`` (model spec but no subcommand) get a useful
    "expected one of these" message instead of ChimeraX's generic
    ``Unknown command: isolde validate #1.2``.
    '''
    raise UserError(
        "'isolde validate' requires a subcommand. Available: "
        + ', '.join(_VALIDATE_SUBCOMMANDS)
        + ". Example: 'isolde validate clashes #1'."
    )


def _full_list_hint(cmd_name, model_spec, total_count):
    '''
    Build a one-line hint pointing the caller at the ``log`` and
    ``saveFile`` keywords. Returned as a leading-space-prefixed string so
    callers can simply ``summary += _full_list_hint(...)`` without
    worrying about extra whitespace, and as the empty string when there
    is nothing to list.
    '''
    if total_count <= 0:
        return ''
    return (
        " Showing summary only; re-run as '{cmd} {spec} log true' to dump"
        " the full table to the ChimeraX log (capped at {cap} rows), or"
        " add 'saveFile <path>' to write the complete list to disk."
    ).format(cmd=cmd_name, spec=model_spec, cap=_LOG_TABLE_ROW_LIMIT)


def _format_table(columns, rows):
    '''
    Return (header_line, separator_line, [data_line, ...]) as monospace
    aligned strings. Columns are left-justified to the widest cell.
    '''
    str_rows = [[str(c) for c in row] for row in rows]
    widths = [len(c) for c in columns]
    for row in str_rows:
        for i, cell in enumerate(row):
            if len(cell) > widths[i]:
                widths[i] = len(cell)
    fmt = '  '.join('{{:<{}}}'.format(w) for w in widths)
    header_line = fmt.format(*columns)
    sep_line = fmt.format(*['-' * w for w in widths])
    data_lines = [fmt.format(*row) for row in str_rows]
    return header_line, sep_line, data_lines


def _dump_table_to_log(logger, header, columns, rows):
    '''
    Log ``header`` followed by a monospace-aligned ``columns`` / ``rows``
    table to the ChimeraX HTML log, matching the pattern used by the
    ChimeraX ``clashes`` and ``hbonds`` commands.
    '''
    header_line, sep_line, data_lines = _format_table(columns, rows)
    truncated = False
    if len(data_lines) > _LOG_TABLE_ROW_LIMIT:
        data_lines = data_lines[:_LOG_TABLE_ROW_LIMIT]
        truncated = True
    body = '\n'.join([header_line, sep_line] + data_lines)
    if truncated:
        body += ('\n... (log output truncated at {} rows; use the saveFile '
            'option to capture the full table).'.format(_LOG_TABLE_ROW_LIMIT))
    logger.info('<pre>' + header + '\n' + body + '</pre>', is_html=True)


def _write_results_file(path, *, summary, columns, rows, json_payload):
    '''
    Write the full validation table to ``path``. ``.json`` paths get a
    structured JSON dump (the full ``json_payload``); any other extension
    gets a plain UTF-8 text table via ``chimerax.io.open_output`` (same
    helper the ChimeraX ``clashes`` / ``hbonds`` commands use).
    '''
    import json as _json
    from chimerax.io import open_output
    if str(path).lower().endswith('.json'):
        with open_output(path, 'utf-8') as f:
            _json.dump(json_payload, f, indent=2, default=str)
        return
    header_line, sep_line, data_lines = _format_table(columns, rows)
    with open_output(path, 'utf-8') as f:
        f.write(summary + '\n\n')
        f.write(header_line + '\n')
        f.write(sep_line + '\n')
        for line in data_lines:
            f.write(line + '\n')


def _maybe_limit(items, limit):
    '''
    Truncate ``items`` to ``limit`` entries for the inline return value.
    Returns ``(returned_items, truncated, total_count, returned_count)``.
    '''
    total = len(items)
    if limit is None or limit < 0 or total <= limit:
        return items, False, total, total
    return items[:limit], True, total, int(limit)


def classify_peptide_bonds(pdm, residues, *,
        cis_cutoff=None, twisted_delta=None):
    '''
    Classify the peptide bonds in ``residues`` as cis or twisted, using
    ``pdm.get_dihedrals(residues, 'omega')``. Returns a list of dicts -
    one per cis or twisted bond, in the order returned by
    ``pdm.get_dihedrals`` - of the form::

        {'omega': ProperDihedral,
         'res1': Residue, 'res2': Residue,
         'omega_deg': float,
         'is_cis': bool, 'is_twisted': bool,
         'is_proline': bool}

    where:

    - ``omega`` is the omega :class:`ProperDihedral` object;
    - ``res1`` / ``res2`` are the N- and C-terminal residues of the bond;
    - ``omega_deg`` is the signed omega angle in degrees;
    - ``is_cis`` is ``|omega| < cis_cutoff``;
    - ``is_twisted`` is ``cis_cutoff <= |omega| < pi - twisted_delta``;
    - ``is_proline`` is true when ``res2.name == 'PRO'``. Cis-Pro bonds are
      valid biology; callers typically split them out from cis non-Pro.

    Cutoffs default to ``isolde.constants.defaults.CIS_PEPTIDE_BOND_CUTOFF``
    and ``defaults.TWISTED_PEPTIDE_BOND_DELTA`` (both expressed in radians).
    The classifier is shared by ``isolde validate peptidebonds`` and
    ISOLDE's GUI "Peptide Bond Validation" panel.
    '''
    import numpy
    from math import pi
    from ..constants import defaults
    if cis_cutoff is None:
        cis_cutoff = defaults.CIS_PEPTIDE_BOND_CUTOFF
    if twisted_delta is None:
        twisted_delta = defaults.TWISTED_PEPTIDE_BOND_DELTA

    omegas = pdm.get_dihedrals(residues, 'omega')
    angles = omegas.angles
    abs_angles = numpy.abs(angles)
    cis_mask = abs_angles < cis_cutoff
    twisted_mask = numpy.logical_and(
        abs_angles >= cis_cutoff, abs_angles < pi - twisted_delta)
    iffy_mask = numpy.logical_or(cis_mask, twisted_mask)
    iffy_omegas = omegas[iffy_mask]
    iffy_angles_deg = numpy.degrees(angles[iffy_mask])
    iffy_cis_mask = cis_mask[iffy_mask]
    iffy_twisted_mask = twisted_mask[iffy_mask]

    out = []
    for omega, angle_deg, is_cis, is_twisted in zip(
            iffy_omegas, iffy_angles_deg, iffy_cis_mask, iffy_twisted_mask):
        res1, res2 = omega.atoms.unique_residues
        out.append({
            'omega': omega,
            'res1': res1,
            'res2': res2,
            'omega_deg': float(angle_deg),
            'is_cis': bool(is_cis),
            'is_twisted': bool(is_twisted),
            'is_proline': res2.name == 'PRO',
        })
    return out


def _compute_peptide_bond_report(session, structure):
    '''
    Compute the structured cis / twisted peptide-bond report for a single
    ``structure``. Pure compute - no logging, no file output, no inline
    truncation. Returns a dict with the same shape as
    ``isolde_validate_peptidebonds`` (minus pagination fields):

        {'model', 'n_residues', 'n_cis_nonpro', 'n_cis_pro',
         'n_twisted', 'n_iffy', 'items'}

    Each item carries ``chain_id``, ``residue_pair_label``, the two
    residue summaries (``res1``, ``res2``), ``omega_deg``,
    ``conformation`` (``'cis'`` or ``'twisted'``) and ``is_proline``.
    Items are sorted with twisted bonds first, then by chain / number.
    '''
    from chimerax.atomic import Residue
    from ..session_extensions import get_proper_dihedral_mgr

    pdm = get_proper_dihedral_mgr(session)
    aa_residues = structure.residues[
        structure.residues.polymer_types == Residue.PT_AMINO]
    raw_items = classify_peptide_bonds(pdm, aa_residues)

    items = []
    n_cis_pro = 0
    n_cis_nonpro = 0
    n_twisted = 0
    for it in raw_items:
        res1 = it['res1']
        res2 = it['res2']
        if it['is_cis']:
            conformation = 'cis'
            if it['is_proline']:
                n_cis_pro += 1
            else:
                n_cis_nonpro += 1
        else:
            conformation = 'twisted'
            n_twisted += 1
        items.append({
            'chain_id': res2.chain_id,
            'residue_pair_label': '{}:{}-{}:{}'.format(
                res1.name, res1.number, res2.name, res2.number),
            'res1': _residue_summary(res1),
            'res2': _residue_summary(res2),
            'omega_deg': float(it['omega_deg']),
            'conformation': conformation,
            'is_proline': bool(it['is_proline']),
        })
    items.sort(key=lambda x: (
        x['conformation'] != 'twisted',
        x['chain_id'],
        x['res2']['number'],
    ))

    return {
        'model': structure.atomspec,
        'n_residues': int(len(aa_residues)),
        'n_cis_nonpro': int(n_cis_nonpro),
        'n_cis_pro': int(n_cis_pro),
        'n_twisted': int(n_twisted),
        'n_iffy': int(len(items)),
        'items': items,
    }


def _compute_rama_report(session, structure, *, include='outliers'):
    '''
    Compute the structured Ramachandran report for a single ``structure``.
    Pure compute - no logging, no file output, no inline truncation.

    ``include`` selects which residues appear in ``items``:
    ``'outliers'`` (default), ``'allowed'`` (outliers + allowed) or
    ``'all'``. Summary counts always cover the full model.
    '''
    import numpy
    from chimerax.atomic import Residue
    from ..session_extensions import get_ramachandran_mgr

    mgr = get_ramachandran_mgr(session)
    RamaBin = mgr.RamaBin
    RamaCase = mgr.RamaCase
    case_names = {
        int(RamaCase.NONE): 'n/a',
        int(RamaCase.CISPRO): 'cis-Pro',
        int(RamaCase.TRANSPRO): 'trans-Pro',
        int(RamaCase.GLYCINE): 'Gly',
        int(RamaCase.PREPRO): 'pre-Pro',
        int(RamaCase.ILEVAL): 'Ile/Val',
        int(RamaCase.GENERAL): 'general',
    }
    bin_names = {
        int(RamaBin.FAVORED): 'favored',
        int(RamaBin.ALLOWED): 'allowed',
        int(RamaBin.OUTLIER): 'outlier',
        int(RamaBin.NA): 'n/a',
    }

    include_choice = (include or 'outliers').lower()
    include_sets = {
        'outliers': {int(RamaBin.OUTLIER)},
        'allowed': {int(RamaBin.OUTLIER), int(RamaBin.ALLOWED)},
        'all': {int(RamaBin.OUTLIER), int(RamaBin.ALLOWED),
            int(RamaBin.FAVORED)},
    }
    if include_choice not in include_sets:
        raise UserError(
            "include must be one of 'outliers', 'allowed', 'all' "
            "(got {!r}).".format(include))
    keep_bins = include_sets[include_choice]

    aa_residues = structure.residues[
        structure.residues.polymer_types == Residue.PT_AMINO]
    ramas = mgr.get_ramas(aa_residues)
    scores, cases = mgr.validate(ramas)
    bins = mgr.bin_scores(scores, cases)
    not_na = bins != int(RamaBin.NA)
    valid_ramas = ramas[not_na]
    scores = scores[not_na]
    cases = cases[not_na]
    bins = bins[not_na]
    phipsis_deg = (numpy.degrees(valid_ramas.phipsis)
        if len(valid_ramas) else numpy.zeros((0, 2)))
    residues = valid_ramas.residues

    n_favored = int((bins == int(RamaBin.FAVORED)).sum())
    n_allowed = int((bins == int(RamaBin.ALLOWED)).sum())
    n_outlier = int((bins == int(RamaBin.OUTLIER)).sum())
    n_total_scored = int(len(valid_ramas))

    severity_rank = {'outlier': 0, 'allowed': 1, 'favored': 2, 'n/a': 3}
    items = []
    for i in range(n_total_scored):
        bin_i = int(bins[i])
        if bin_i not in keep_bins:
            continue
        r = residues[i]
        phi, psi = phipsis_deg[i]
        items.append({
            'chain_id': r.chain_id,
            'name': r.name,
            'number': int(r.number),
            'insertion_code': r.insertion_code,
            'spec': r.atomspec,
            'phi_deg': float(phi),
            'psi_deg': float(psi),
            'score': float(scores[i]),
            'classification': bin_names[bin_i],
            'case': case_names[int(cases[i])],
        })
    items.sort(key=lambda it: (
        severity_rank.get(it['classification'], 9),
        it['chain_id'],
        it['number'],
    ))

    return {
        'model': structure.atomspec,
        'include': include_choice,
        'n_scorable': n_total_scored,
        'n_favored': n_favored,
        'n_allowed': n_allowed,
        'n_outlier': n_outlier,
        'items': items,
    }


def _compute_rotamer_report(session, structure, *, include='nonfavored'):
    '''
    Compute the structured rotamer report for a single ``structure``.
    Pure compute - no logging, no file output, no inline truncation.

    ``include`` selects which residues appear in ``items``: ``'outliers'``,
    ``'nonfavored'`` (default; outliers + allowed) or ``'all'``. Summary
    counts always cover the full set of rotameric residues.
    '''
    import numpy
    from ..session_extensions import get_rotamer_mgr

    mgr = get_rotamer_mgr(session)
    allowed_cutoff, outlier_cutoff = mgr.cutoffs

    include_choice = (include or 'nonfavored').lower()
    include_sets = {
        'outliers': {'outlier'},
        'nonfavored': {'outlier', 'allowed'},
        'all': {'outlier', 'allowed', 'favored'},
    }
    if include_choice not in include_sets:
        raise UserError(
            "include must be one of 'outliers', 'nonfavored', 'all' "
            "(got {!r}).".format(include))
    keep = include_sets[include_choice]

    rotas = mgr.get_rotamers(structure.residues)
    n_total = len(rotas)
    if n_total:
        scores = mgr.validate_rotamers(rotas)
        residues = rotas.residues
    else:
        scores = numpy.zeros(0, dtype=float)
        from chimerax.atomic import Residues
        residues = Residues()

    favored_mask = scores >= allowed_cutoff
    outlier_mask = scores < outlier_cutoff
    allowed_mask = numpy.logical_and(~favored_mask, ~outlier_mask)
    n_favored = int(favored_mask.sum())
    n_allowed = int(allowed_mask.sum())
    n_outlier = int(outlier_mask.sum())

    severity_rank = {'outlier': 0, 'allowed': 1, 'favored': 2}
    items = []
    for i in range(n_total):
        s = float(scores[i])
        if s >= allowed_cutoff:
            cls = 'favored'
        elif s >= outlier_cutoff:
            cls = 'allowed'
        else:
            cls = 'outlier'
        if cls not in keep:
            continue
        r = residues[i]
        items.append({
            'chain_id': r.chain_id,
            'name': r.name,
            'number': int(r.number),
            'insertion_code': r.insertion_code,
            'spec': r.atomspec,
            'resname': r.name,
            'score': s,
            'classification': cls,
        })
    items.sort(key=lambda it: (
        severity_rank[it['classification']],
        it['chain_id'],
        it['number'],
    ))

    return {
        'model': structure.atomspec,
        'include': include_choice,
        'cutoff_allowed': float(allowed_cutoff),
        'cutoff_outlier': float(outlier_cutoff),
        'n_rotameric': int(n_total),
        'n_favored': n_favored,
        'n_allowed': n_allowed,
        'n_outlier': n_outlier,
        'items': items,
    }


def isolde_validate_peptidebonds(session, model=None,
        save_file=None, log=False, limit=None):
    '''
    Report cis and twisted peptide bonds in ``model`` (or ISOLDE's currently
    selected model), using the same omega-dihedral classification that
    ISOLDE's "Peptide Bond Validation" panel applies. A bond is classified
    ``cis`` when |omega| < CIS_PEPTIDE_BOND_CUTOFF (default 30 deg) and
    ``twisted`` when |omega| is between that cutoff and
    pi - TWISTED_PEPTIDE_BOND_DELTA. Cis-proline bonds are valid and are
    reported separately from cis non-proline bonds.

    Read only - the model is never modified.

    Returns
    -------
    dict
        Summary counts (``n_residues``, ``n_cis_nonpro``, ``n_cis_pro``,
        ``n_twisted``, ``n_iffy``) plus a per-bond ``items`` list. Each
        item carries the chain, both residues, omega angle in degrees,
        ``conformation`` (``cis`` / ``twisted``), and an ``is_proline``
        flag for the C-terminal residue. ``truncated`` indicates whether
        ``limit`` clipped the list (``total_count`` is the unclipped size).
    '''
    m = _resolve_model(session, model)
    logger = session.logger

    data = _compute_peptide_bond_report(session, m)
    items = data['items']

    columns = ['Chain', 'Residues', 'Omega (deg)', 'Conformation', 'Is Pro?']
    rows = [
        (
            it['chain_id'],
            it['residue_pair_label'],
            '{:.1f}'.format(it['omega_deg']),
            it['conformation'],
            'yes' if it['is_proline'] else 'no',
        ) for it in items
    ]

    summary = (
        'ISOLDE peptide bond check ({}): {} of {} amino-acid residues '
        'are cis or twisted ({} cis non-Pro, {} cis Pro, {} twisted).'.format(
            m.atomspec, data['n_iffy'], data['n_residues'],
            data['n_cis_nonpro'], data['n_cis_pro'], data['n_twisted'],
        )
    )
    summary += _full_list_hint(
        'isolde validate peptidebonds', m.atomspec, data['n_iffy'])
    if data['n_cis_nonpro'] > 0 or data['n_twisted'] > 0:
        logger.warning(summary)
    else:
        logger.info(summary)

    if log:
        _dump_table_to_log(logger, summary, columns, rows)

    returned_items, truncated, total_count, returned_count = _maybe_limit(items, limit)

    result = dict(data)
    result['items'] = returned_items
    result['returned_count'] = int(returned_count)
    result['total_count'] = int(total_count)
    result['truncated'] = bool(truncated)

    if save_file is not None:
        _write_results_file(save_file,
            summary=summary, columns=columns, rows=rows,
            json_payload=dict(data, items=items))
        logger.info('Wrote peptide bond report to {}'.format(save_file))

    return result


def isolde_validate_rama(session, model=None, include='outliers',
        save_file=None, log=False, limit=None):
    '''
    Report Ramachandran scoring for protein residues in ``model`` (or
    ISOLDE's currently selected model), using the same MolProbity contours
    and bin cutoffs that ISOLDE's Ramachandran plot and validator use.

    ``include`` selects which residues appear in the per-residue list:
    ``'outliers'`` (default), ``'allowed'`` (outliers + allowed) or
    ``'all'`` (favored too). Summary counts always cover the full model.

    Read only - the model is never modified, and no live annotators are
    created. To toggle the live 3D annotators see the existing ``rama``
    command.

    Returns
    -------
    dict
        Summary counts (``n_scorable``, ``n_favored``, ``n_allowed``,
        ``n_outlier``) plus a per-residue ``items`` list with phi/psi in
        degrees, the MolProbity ``score``, ``classification`` and the
        Ramachandran ``case`` (e.g. ``general``, ``Gly``, ``trans-Pro``).
    '''
    m = _resolve_model(session, model)
    logger = session.logger

    data = _compute_rama_report(session, m, include=include)
    items = data['items']

    columns = ['Chain', 'Residue', 'Phi', 'Psi', 'Score', 'Class', 'Case']
    rows = [
        (
            it['chain_id'],
            '{} {}'.format(it['name'], it['number']),
            '{:.1f}'.format(it['phi_deg']),
            '{:.1f}'.format(it['psi_deg']),
            '{:.4f}'.format(it['score']),
            it['classification'],
            it['case'],
        ) for it in items
    ]

    summary = (
        'ISOLDE Ramachandran check ({}): {} outliers, {} allowed, {} favored '
        '(of {} scorable residues).'.format(
            m.atomspec, data['n_outlier'], data['n_allowed'],
            data['n_favored'], data['n_scorable'],
        )
    )
    summary += _full_list_hint(
        'isolde validate rama', m.atomspec, len(items))
    if data['n_outlier'] > 0:
        logger.warning(summary)
    else:
        logger.info(summary)

    if log:
        _dump_table_to_log(logger, summary, columns, rows)

    returned_items, truncated, total_count, returned_count = _maybe_limit(items, limit)

    result = dict(data)
    result['items'] = returned_items
    result['returned_count'] = int(returned_count)
    result['total_count'] = int(total_count)
    result['truncated'] = bool(truncated)

    if save_file is not None:
        _write_results_file(save_file,
            summary=summary, columns=columns, rows=rows,
            json_payload=dict(data, items=items))
        logger.info('Wrote Ramachandran report to {}'.format(save_file))

    return result


def isolde_validate_rotamers(session, model=None, include='nonfavored',
        save_file=None, log=False, limit=None):
    '''
    Report rotamer scoring for sidechain-bearing residues in ``model`` (or
    ISOLDE's currently selected model), using the same MolProbity contours
    that ISOLDE's "Rotamer Validation" panel applies. Each residue is
    classified ``outlier`` / ``allowed`` / ``favored`` against the current
    P-value cutoffs on the session ``RotaMgr``.

    ``include`` selects which residues appear in the per-residue list:
    ``'nonfavored'`` (default; outliers + allowed), ``'outliers'`` or
    ``'all'``. Summary counts always cover the full set of rotameric
    residues.

    Read only - the model is never modified, and no live annotators are
    created. To toggle the live 3D annotators see the existing ``rota``
    command.

    Returns
    -------
    dict
        Summary counts (``n_rotameric``, ``n_favored``, ``n_allowed``,
        ``n_outlier``), cutoff values, and a per-residue ``items`` list
        with the P-value ``score`` and ``classification``.
    '''
    m = _resolve_model(session, model)
    logger = session.logger

    data = _compute_rotamer_report(session, m, include=include)
    items = data['items']

    columns = ['Chain', 'Residue', 'Resname', 'P', 'Class']
    rows = [
        (
            it['chain_id'],
            str(it['number']),
            it['resname'],
            '{:.4f}'.format(it['score']),
            it['classification'],
        ) for it in items
    ]

    summary = (
        'ISOLDE rotamer check ({}): {} outliers, {} allowed, {} favored '
        '(of {} rotameric residues; cutoffs allowed>={:.3f}, '
        'outlier<{:.3f}).'.format(
            m.atomspec, data['n_outlier'], data['n_allowed'],
            data['n_favored'], data['n_rotameric'],
            data['cutoff_allowed'], data['cutoff_outlier'],
        )
    )
    summary += _full_list_hint(
        'isolde validate rotamers', m.atomspec, len(items))
    if data['n_outlier'] > 0:
        logger.warning(summary)
    else:
        logger.info(summary)

    if log:
        _dump_table_to_log(logger, summary, columns, rows)

    returned_items, truncated, total_count, returned_count = _maybe_limit(items, limit)

    result = dict(data)
    result['items'] = returned_items
    result['returned_count'] = int(returned_count)
    result['total_count'] = int(total_count)
    result['truncated'] = bool(truncated)

    if save_file is not None:
        _write_results_file(save_file,
            summary=summary, columns=columns, rows=rows,
            json_payload=dict(data, items=items))
        logger.info('Wrote rotamer report to {}'.format(save_file))

    return result


def isolde_validate_clashes(session, model=None,
        save_file=None, log=False, limit=200):
    '''
    Report steric clashes in ``model`` (or ISOLDE's currently selected
    model), using ISOLDE's ``unique_clashes`` wrapper around the ChimeraX
    ``clashes`` machinery. Each clash carries the two atoms, the van der
    Waals overlap in Angstroms, and a ``severity`` of either ``strict``
    (overlap >= STRICT_CUTOFF) or ``severe`` (overlap >= SEVERE_CUTOFF).

    Read only - the model is never modified. Returns a dict with summary
    counts and a per-clash ``items`` list, defaulting to the worst 200
    clashes inline (use ``limit`` to widen, or ``save_file`` to capture
    the full list).

    Returns
    -------
    dict
        Summary counts (``n_total``, ``n_severe``, ``n_strict``), cutoff
        values, and a per-clash ``items`` list sorted by descending
        overlap.
    '''
    m = _resolve_model(session, model)
    logger = session.logger

    from .clashes import (
        unique_clashes, clash_atom_label, STRICT_CUTOFF, SEVERE_CUTOFF,
    )
    clashes = unique_clashes(session, m.atoms)
    n_total = int(len(clashes))

    items = []
    n_severe = 0
    n_strict = 0
    for clash in clashes:
        a1, a2 = clash.atoms
        overlap = float(clash.overlap)
        if overlap >= SEVERE_CUTOFF:
            severity = 'severe'
            n_severe += 1
        else:
            severity = 'strict'
            n_strict += 1
        items.append({
            'atom1_spec': a1.atomspec,
            'atom1_label': clash_atom_label(a1),
            'atom2_spec': a2.atomspec,
            'atom2_label': clash_atom_label(a2),
            'overlap': overlap,
            'severity': severity,
        })

    columns = ['Atom 1', 'Atom 2', 'Overlap', 'Severity']
    rows = [
        (
            it['atom1_label'],
            it['atom2_label'],
            '{:.2f}'.format(it['overlap']),
            it['severity'],
        ) for it in items
    ]

    summary = (
        'ISOLDE clash check ({}): {} unique clashes ({} severe, {} strict; '
        'cutoffs severe>={:.2f}, strict>={:.2f} A).'.format(
            m.atomspec, n_total, n_severe, n_strict,
            float(SEVERE_CUTOFF), float(STRICT_CUTOFF),
        )
    )
    summary += _full_list_hint(
        'isolde validate clashes', m.atomspec, n_total)
    if n_severe > 0:
        logger.warning(summary)
    else:
        logger.info(summary)

    if log:
        _dump_table_to_log(logger, summary, columns, rows)

    returned_items, truncated, total_count, returned_count = _maybe_limit(items, limit)

    result = {
        'model': m.atomspec,
        'severe_cutoff': float(SEVERE_CUTOFF),
        'strict_cutoff': float(STRICT_CUTOFF),
        'n_total': n_total,
        'n_severe': int(n_severe),
        'n_strict': int(n_strict),
        'items': returned_items,
        'returned_count': int(returned_count),
        'total_count': int(total_count),
        'truncated': bool(truncated),
    }

    if save_file is not None:
        _write_results_file(save_file,
            summary=summary, columns=columns, rows=rows,
            json_payload=dict(result, items=items))
        logger.info('Wrote clash report to {}'.format(save_file))

    return result


def register_validate_commands(logger):
    from chimerax.core.commands import (
        register, CmdDesc, BoolArg, IntArg, EnumOf, SaveFileNameArg,
    )
    from ..cmd.argspec import IsoldeStructureArg

    common_kw = [
        ('save_file', SaveFileNameArg),
        ('log', BoolArg),
        ('limit', IntArg),
    ]

    desc_pep = CmdDesc(
        optional=[('model', IsoldeStructureArg)],
        keyword=list(common_kw),
        synopsis='Validation: report cis and twisted peptide bonds.',
    )
    register('isolde validate peptidebonds', desc_pep,
        isolde_validate_peptidebonds, logger=logger)

    desc_rama = CmdDesc(
        optional=[('model', IsoldeStructureArg)],
        keyword=list(common_kw) + [
            ('include', EnumOf(('outliers', 'allowed', 'all'))),
        ],
        synopsis='Validation: report Ramachandran scoring with phi/psi/score.',
    )
    register('isolde validate rama', desc_rama,
        isolde_validate_rama, logger=logger)

    desc_rota = CmdDesc(
        optional=[('model', IsoldeStructureArg)],
        keyword=list(common_kw) + [
            ('include', EnumOf(('outliers', 'nonfavored', 'all'))),
        ],
        synopsis='Validation: report rotamer scoring (favored/allowed/outlier).',
    )
    register('isolde validate rotamers', desc_rota,
        isolde_validate_rotamers, logger=logger)

    desc_clash = CmdDesc(
        optional=[('model', IsoldeStructureArg)],
        keyword=list(common_kw),
        synopsis='Validation: report steric clashes.',
    )
    register('isolde validate clashes', desc_clash,
        isolde_validate_clashes, logger=logger)

    # Parent command: catches ``isolde validate`` (no subcommand) and
    # ``isolde validate <model>`` (model spec but no subcommand) and
    # turns them into a helpful "expected one of: ..." error instead of
    # ChimeraX's generic ``Unknown command``. The ``optional`` model arg
    # exists so the parser cleanly consumes the spec - the handler always
    # raises before doing any work.
    desc_top = CmdDesc(
        optional=[('model', IsoldeStructureArg)],
        synopsis='Run an ISOLDE validator (requires a subcommand: {}).'.format(
            ', '.join(_VALIDATE_SUBCOMMANDS)),
    )
    register('isolde validate', desc_top, isolde_validate, logger=logger)
