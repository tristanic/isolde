# @Author: Tristan Croll
# @Date:   27-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 27-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Unit check for the human-facing announcement path (REST/MCP transparency).

Run inside ChimeraX (no GUI / no simulation needed):

    run_chimerax.bat --nogui --exit --script src/tests/test_announce.py

A remote agent's commands run inside a StringPlainTextLog whose
``excludes_other_logs`` is True, so ordinary log output is captured for the
remote and HIDDEN from every other log (in the GUI, that means the human sees
nothing). ``announce_to_user`` must punch through that exclusion for one message
so a human sharing the session can follow the agent -- while the capture keeps a
copy and is left intact afterwards. This reproduces that exact log topology with
a recording log standing in for the GUI log, and asserts:

  - under an active exclusive capture, a normal logger.info is EXCLUDED from the
    recorder (the bug we are guarding against);
  - announce_to_user reaches BOTH the recorder and the capture;
  - the capture's excludes_other_logs flag is restored to True afterwards;
  - _outcome_summary renders each witness shape into a short human phrase, and
    stays silent ('') when there is nothing worth announcing.
'''
from chimerax.core.logger import PlainTextLog, StringPlainTextLog


def _fail(msg):
    print('FAIL:', msg)
    raise SystemExit(1)


class _Recorder(PlainTextLog):
    '''Stand-in for the GUI log: does NOT exclude other logs, records messages.'''
    excludes_other_logs = False

    def __init__(self):
        self.messages = []

    def log(self, level, msg):
        self.messages.append(msg)
        return True   # 'displayed'

    def text(self):
        return ''.join(self.messages)


def run(session):
    from chimerax.isolde.cmd.invoke import announce_to_user, _outcome_summary

    logger = session.logger
    recorder = _Recorder()
    logger.add_log(recorder)
    try:
        with StringPlainTextLog(logger) as cap:
            # Control: an ordinary info() is consumed by the exclusive capture
            # and never reaches the recorder -- the exact "human sees nothing"
            # condition we are fixing.
            logger.info('CONTROL-LINE')
            if 'CONTROL-LINE' not in cap.getvalue():
                _fail('capture did not receive the control line')
            if 'CONTROL-LINE' in recorder.text():
                _fail('control line leaked to recorder -- capture not exclusive?')
            print('PASS: under exclusive capture, plain info() is hidden from the GUI log')

            # announce_to_user must reach the recorder (GUI) too...
            announce_to_user(session, 'ANNOUNCE-LINE')
            if 'ANNOUNCE-LINE' not in recorder.text():
                _fail('announce_to_user did NOT reach the GUI/recorder log')
            print('PASS: announce_to_user reaches the GUI log despite the capture')

            # ...and the agent's captured remote log gets a copy.
            if 'ANNOUNCE-LINE' not in cap.getvalue():
                _fail('announce_to_user did not also land in the captured remote log')
            print('PASS: the captured remote log keeps a copy of the announcement')

            # ...and the exclusion must be restored so subsequent command output
            # is captured normally again.
            if not cap.excludes_other_logs:
                _fail('announce_to_user left excludes_other_logs cleared')
            logger.info('AFTER-LINE')
            if 'AFTER-LINE' in recorder.text():
                _fail('capture exclusivity not restored (AFTER-LINE leaked)')
            print('PASS: capture exclusivity restored after the announcement')
    finally:
        logger.remove_log(recorder)

    # _outcome_summary: each witness shape -> a short human phrase.
    cases = [
        ({'simulation_running': True, 'cmd': 'start', 'started_ok': True}, 'simulation started'),
        ({'simulation_running': False, 'cmd': 'stop', 'stopped': True}, 'simulation stopped'),
        ({'simulation_running': True, 'cmd': 'stop', 'stopped': False}, 'simulation stopping'),
        ({'restraint_total_delta': 12, 'restraint_enabled_delta': 12}, '+12 restraints'),
        ({'restraint_total_delta': 0, 'restraint_enabled_delta': -8}, '8 restraints released'),
        ({'launched': True}, 'refinement launched (runs in the background)'),
        ({'rmsd_moved': 0.43, 'max_atom_shift': 1.20}, 'moved RMSD 0.43'),
        ({'rmsd_moved': 0.0, 'sim_running': True}, 'target set'),
    ]
    for witness, expect in cases:
        got = _outcome_summary(witness)
        if expect not in got:
            _fail('outcome summary for %s was %r, expected to contain %r'
                  % (witness, got, expect))
    # Silent when nothing happened / read-only.
    if _outcome_summary({}) != '':
        _fail('empty witness should summarise to ""')
    if _outcome_summary({'rmsd_moved': 0.0}) != 'no movement':
        _fail('zero-move (no sim) should be "no movement"')
    print('PASS: _outcome_summary renders every witness shape (and stays silent on none)')

    print('ALL PASS')


try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
