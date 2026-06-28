# @Author: Tristan Croll
# @Date:   28-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 28-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Unit check for the ISOLDE philosophy layer (agent/MCP control surface).

Run inside ChimeraX (no GUI / no simulation needed):

    run_chimerax.bat --nogui --exit --script src/tests/test_philosophy.py

Asserts:
  - the canonical module is well-formed (10 principles, unique 1..10 ids) and the
    condensed INSTRUCTIONS quote every principle title (so they can't drift);
  - each soft-guardrail detector fires on the violating shape and STAYS SILENT on
    the benign one (force-fit only on a coupling *increase*; strong-tug only ABOVE
    default; large-move only past the threshold; clean-metric only at zero outliers);
  - announce_high_severity punches a high-severity flag through to the GUI log even
    under an exclusive REST capture (the human-in-the-loop safety net, principle 7),
    and emits nothing for warn/info flags;
  - every surface reads from the one module: the agent_philosophy REST method
    returns the module's payload, and agent_tools() carries MANIFEST_PREAMBLE.
'''
from chimerax.core.logger import PlainTextLog, StringPlainTextLog


def _fail(msg):
    print('FAIL:', msg)
    raise SystemExit(1)


def _eq(got, want, what):
    if got != want:
        _fail('%s: got %r, expected %r' % (what, got, want))


class _Recorder(PlainTextLog):
    '''Stand-in for the GUI log: records messages, does not exclude other logs.'''
    excludes_other_logs = False

    def __init__(self):
        self.messages = []

    def log(self, level, msg):
        self.messages.append(msg)
        return True

    def text(self):
        return ''.join(self.messages)


def run(session):
    from chimerax.isolde.cmd import agent_philosophy as ap

    # --- canonical content -------------------------------------------------
    ids = [p['id'] for p in ap.PRINCIPLES]
    _eq(len(ap.PRINCIPLES), 10, 'number of principles')
    _eq(ids, list(range(1, 11)), 'principle ids 1..10 in order')
    for p in ap.PRINCIPLES:
        if p['title'] not in ap.INSTRUCTIONS:
            _fail('INSTRUCTIONS missing principle title: %r' % p['title'])
        if p['title'] not in ap.principles_document() or p['text'] not in ap.principles_document():
            _fail('principles_document missing principle %d text/title' % p['id'])
    print('PASS: 10 principles; INSTRUCTIONS + document quote every title/text')

    # --- operational layer: workflow + rules of thumb ----------------------
    if not ap.WORKFLOW or not ap.RULES_OF_THUMB:
        _fail('WORKFLOW / RULES_OF_THUMB must be non-empty')
    doc = ap.principles_document()
    for needle in ('Recommended workflow', 'Rules of thumb', 'problem_zones',
                   'pLDDT', 'register-shifter'):
        if needle not in doc:
            _fail('principles_document missing operational content: %r' % needle)
    # The operational depth belongs in the on-demand resource, NOT the always-in-
    # context constitution -- keep INSTRUCTIONS tight.
    if 'Recommended workflow' in ap.INSTRUCTIONS or 'pLDDT' in ap.INSTRUCTIONS:
        _fail('INSTRUCTIONS should stay the tight constitution (no workflow/rules inlined)')
    print('PASS: workflow + rules of thumb in the document, kept out of INSTRUCTIONS')

    # --- guardrail detectors: fire on violation, silent on benign ----------
    def one(flags, principle, severity, what):
        if len(flags) != 1:
            _fail('%s: expected exactly one flag, got %r' % (what, flags))
        _eq(flags[0]['principle'], principle, what + ' principle')
        _eq(flags[0]['severity'], severity, what + ' severity')

    one(ap.flag_force_fit(1.0, 2.0), 4, 'high', 'force-fit (increase)')
    _eq(ap.flag_force_fit(2.0, 2.0), [], 'force-fit (unchanged) silent')
    _eq(ap.flag_force_fit(2.0, 1.0), [], 'force-fit (decrease) silent')
    _eq(ap.flag_force_fit(None, 5.0), [], 'force-fit (no prior) silent')

    one(ap.flag_strong_tug(1000.0, 500.0), 4, 'warn', 'strong-tug (above default)')
    _eq(ap.flag_strong_tug(500.0, 500.0), [], 'strong-tug (default) silent')
    _eq(ap.flag_strong_tug(100.0, 500.0), [], 'strong-tug (weaker) silent')

    big = ap.LARGE_MOVE_ANGSTROM + 1.0
    one(ap.flag_large_move({'max_atom_shift': big}), 5, 'warn', 'large-move')
    _eq(ap.flag_large_move({'max_atom_shift': 0.2}), [], 'small-move silent')
    _eq(ap.flag_large_move({'rmsd_moved': 0.1}), [], 'no-shift-key silent')
    _eq(ap.flag_large_move(None), [], 'no-witness silent')

    one(ap.flag_clean_validation({'n_outlier': 0}), 2, 'info', 'clean-metric (rama)')
    one(ap.flag_clean_validation({'n_twisted': 0, 'n_iffy': 2}), 2, 'info', 'clean-metric (peptide)')
    _eq(ap.flag_clean_validation({'n_outlier': 3}), [], 'dirty-metric silent')
    _eq(ap.flag_clean_validation({'count': 5}), [], 'no-count-key silent')

    # evaluate_guardrails (the invoke choke-point entry)
    one(ap.evaluate_guardrails(category='validation', result={'n_outlier': 0}),
        2, 'info', 'evaluate: clean validation')
    one(ap.evaluate_guardrails(category='manipulation',
                               witness={'max_atom_shift': big}),
        5, 'warn', 'evaluate: large move')
    _eq(ap.evaluate_guardrails(category='query', result={'n_outlier': 0}), [],
        'evaluate: clean count under non-validation category is silent')
    _eq(ap.evaluate_guardrails(category='manipulation',
                               witness={'max_atom_shift': 0.1}), [],
        'evaluate: small move silent')
    print('PASS: every guardrail fires on violation and stays silent on benign')

    # --- announce_high_severity reaches the GUI log under exclusive capture -
    logger = session.logger
    recorder = _Recorder()
    logger.add_log(recorder)
    try:
        with StringPlainTextLog(logger) as cap:
            ap.announce_high_severity(session, ap.flag_force_fit(1.0, 2.0))
            if 'force-fit' not in recorder.text().lower():
                _fail('high-severity flag did not reach the GUI/recorder log')
            if 'force-fit' not in cap.getvalue().lower():
                _fail('high-severity announce did not also land in the remote capture')
            before = recorder.text()
            ap.announce_high_severity(session, ap.flag_strong_tug(1000.0, 500.0))
            ap.announce_high_severity(session, [])
            if recorder.text() != before:
                _fail('warn/empty flags should NOT be announced to the GUI log')
    finally:
        logger.remove_log(recorder)
    print('PASS: high-severity flags reach the GUI log; warn/info stay quiet')

    # --- every surface reads from the one module ---------------------------
    from chimerax.isolde.remote_control.rest_server.agent_methods import (
        agent_philosophy, agent_tools)
    payload = agent_philosophy(session)
    _eq(payload['instructions'], ap.INSTRUCTIONS, 'REST agent_philosophy instructions')
    _eq(payload['preamble'], ap.MANIFEST_PREAMBLE, 'REST agent_philosophy preamble')
    _eq(len(payload['principles']), 10, 'REST agent_philosophy principle count')
    _eq(payload['workflow'], list(ap.WORKFLOW), 'REST agent_philosophy workflow')
    _eq(payload['rules_of_thumb'], list(ap.RULES_OF_THUMB), 'REST agent_philosophy rules')
    _eq(agent_tools(session).get('preamble'), ap.MANIFEST_PREAMBLE,
        'agent_tools manifest preamble')
    print('PASS: REST agent_philosophy + agent_tools read from the canonical module')

    print('ALL PASS')


try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
