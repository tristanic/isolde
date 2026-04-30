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
        from chimerax.atomic import Residues, concatenate
        residues = concatenate([m.residues for m in structures])
        mgr = sx.get_rotamer_mgr(session)
        rotamers = mgr.get_rotamers(residues)
        report_str = 'NON-FAVOURED ROTAMERS: \n'
        nf, scores = mgr.non_favored_rotamers(rotamers)
        for r, score in zip(nf, scores):
            report_str += '#{:<6} {}:\t{} {} (P={:.4f})\n'.format(
                r.residue.structure.id_string, r.residue.chain_id, r.residue.name,
                r.residue.number, score
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
        from chimerax.atomic import Residues, concatenate
        residues = concatenate([m.residues for m in structures])
        mgr = sx.get_ramachandran_mgr(session)
        report_str = 'RAMACHANDRAN OUTLIERS: \n'
        outliers = mgr.outliers(residues)
        for outlier in outliers:
            report_str +='#{:<6} {}:\t{} {}\n'.format(
                outlier.structure.id_string, outlier.chain_id, outlier.name, outlier.number)
        report_str += '\nCIS PEPTIDE BONDS: \n'
        cispeps = mgr.cis(residues)
        for cis in cispeps:
            report_str +='#{:<6} {}:\t{} {}\n'.format(
                cis.structure.id_string, cis.chain_id, cis.name, cis.number
            )
        report_str += '\nTWISTED PEPTIDE BONDS: \n'
        twisteds = mgr.twisted(residues)
        for twisted, angle in twisteds:
            report_str += '#{:<6} {}:\t{} {} ({:.1f}°)\n'.format(
                twisted.structure.id_string, twisted.chain_id, twisted.name,
                twisted.number, angle
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
# Read-only ``isolde validate`` commands.
#
# These commands let an agent (or scripted user) ask whether a model is ready
# for an ISOLDE simulation without actually trying to start one. They are safe
# to call at any time: they never construct an OpenMM ``Context``, never raise
# on missing parameters, and never modify the model.
# ---------------------------------------------------------------------------

from chimerax.core.errors import UserError

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


def isolde_validate_hydrogens(session, model=None):
    '''
    Check whether ``model`` (or ISOLDE's currently selected model) appears to
    have a complete set of hydrogens, using the same heuristics ISOLDE's
    "Unparametrised residues" panel applies before launching a simulation.

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


def isolde_validate_parameters(session, model=None, forcefield=None,
        ignore_external_bonds=True):
    '''
    Run ISOLDE's MD-template assignment over the given model (or ISOLDE's
    currently selected model) without starting a simulation, and report any
    residues that are unparametrised or ambiguous.

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


def register_validate_commands(logger):
    from chimerax.core.commands import (
        register, CmdDesc, BoolArg, StringArg,
    )
    from ..cmd.argspec import IsoldeStructureArg

    desc_h = CmdDesc(
        optional=[('model', IsoldeStructureArg)],
        synopsis='Check whether the model has plausible hydrogens',
    )
    register('isolde validate hydrogens', desc_h,
        isolde_validate_hydrogens, logger=logger)

    desc_p = CmdDesc(
        optional=[('model', IsoldeStructureArg)],
        keyword=[
            ('forcefield', StringArg),
            ('ignore_external_bonds', BoolArg),
        ],
        synopsis='Check whether all residues match an MD forcefield template',
    )
    register('isolde validate parameters', desc_p,
        isolde_validate_parameters, logger=logger)
