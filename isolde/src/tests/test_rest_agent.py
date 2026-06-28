# @Author: Tristan Croll
# @Date:   27-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 27-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Phase-B smoke test for the revamped REST substrate. Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_rest_agent.py

Validates (no GUI / no simulation needed):
  - bearer-token auth: unauthenticated request rejected (401), authenticated ok;
  - GET manifest lists the agent endpoints and does NOT expose gated 'run';
  - agent_tools returns the auto-generated typed tool manifest;
  - agent_invoke runs a command and returns a structured result;
  - list_models / resolve_spec transport helpers work;
  - render degrades gracefully when headless (no GL context).

Runs the HTTP client in the same process: safe because ChimeraX's nogui event
loop has not started under --script, so thread_safe executes directly in the
server thread (no main-thread marshaling -> no deadlock).
'''
import json
import time
import http.client


def _fail(msg):
    print('FAIL:', msg)
    raise SystemExit(1)


def _request(port, method, body=None, token=None):
    conn = http.client.HTTPConnection('localhost', port, timeout=30)
    headers = {'Content-type': 'application/json'}
    if token is not None:
        headers['Authorization'] = 'Bearer ' + token
    payload = json.dumps(body).encode('utf-8') if body is not None else None
    conn.request(method, '', payload, headers)
    resp = conn.getresponse()
    status = resp.status
    data = resp.read().decode('utf-8')
    conn.close()
    try:
        parsed = json.loads(data)
    except Exception:
        parsed = data
    return status, parsed


def _post(port, token, cmd, **kwargs):
    return _request(port, 'POST', {'cmd': cmd, 'args': [], 'kwargs': kwargs}, token=token)


def run(session):
    from chimerax.isolde.remote_control.rest_server import cmd as rest_cmd

    rest_cmd.start_server(session)
    # Use the raw module global rather than _get_server(), which nulls the
    # server if httpd isn't bound yet (a startup race with the server thread).
    server = rest_cmd._server
    if server is None:
        _fail('server object was not created')

    # Wait for the HTTP server thread to bind its socket.
    for _ in range(200):
        if server.httpd is not None and server.server_address is not None:
            break
        time.sleep(0.05)
    if server.httpd is None:
        _fail('server httpd never bound')
    port = server.port
    token = server.auth_token
    print('server on port', port, 'token len', len(token or ''))
    if not token:
        _fail('no auth token minted')

    # 1. Auth: no token -> 401
    status, _ = _request(port, 'GET', token=None)
    if status != 401:
        _fail('unauthenticated GET should be 401, got %s' % status)
    print('PASS: unauthenticated request rejected (401)')

    # 2. Auth: wrong token -> 401
    status, _ = _request(port, 'GET', token='wrong-token')
    if status != 401:
        _fail('bad-token GET should be 401, got %s' % status)
    print('PASS: bad token rejected (401)')

    # 3. GET manifest with token -> 200, lists agent endpoints, hides gated run
    status, methods = _request(port, 'GET', token=token)
    if status != 200:
        _fail('authenticated GET should be 200, got %s' % status)
    for needed in ('agent_tools', 'agent_invoke', 'list_models', 'resolve_spec'):
        if needed not in methods:
            _fail('manifest missing endpoint %r' % needed)
    if 'run' in methods:
        _fail("'run' should be gated off (allow_run defaulted false)")
    print('PASS: GET manifest authenticated; agent endpoints present; run gated off')

    # 4. agent_tools -> typed tool manifest
    status, result = _post(port, token, 'agent_tools')
    if status != 200 or 'tools' not in result:
        _fail('agent_tools failed: %s %s' % (status, result))
    tools = result['tools']
    names = {t.get('command') for t in tools}
    if 'isolde sim' not in names or 'isolde validate ramachandran' not in names:
        _fail('agent_tools manifest missing expected commands: %s' % sorted(names))
    print('PASS: agent_tools returned %d typed tools' % result['count'])

    # 5. agent_invoke 'isolde status' -> structured result
    status, result = _post(port, token, 'agent_invoke', name='isolde status', args={})
    if status != 200:
        _fail('agent_invoke failed: %s %s' % (status, result))
    if 'error' in result:
        _fail('agent_invoke returned error: %s' % result)
    inner = result.get('result', {})
    if not isinstance(inner, dict) or 'isolde_started' not in inner:
        _fail('agent_invoke result not the status dict: %s' % result)
    print('PASS: agent_invoke isolde status ->', inner)

    # 6. agent_invoke rejects a non-agent_safe command
    status, result = _post(port, token, 'agent_invoke', name='isolde tutorial', args={})
    if 'error' not in result:
        _fail('invoking non-agent_safe command should error, got %s' % result)
    print('PASS: agent_invoke refuses non-agent_safe command (isolde tutorial)')

    # 7. list_models
    status, result = _post(port, token, 'list_models')
    if status != 200 or 'models' not in result:
        _fail('list_models failed: %s %s' % (status, result))
    print('PASS: list_models -> %d models' % result['count'])

    # 8. resolve_spec (no model loaded; should still return a dict, not crash)
    status, result = _post(port, token, 'resolve_spec', spec='#1')
    if status != 200 or 'spec' not in result:
        _fail('resolve_spec failed: %s %s' % (status, result))
    print('PASS: resolve_spec -> valid=%s n_atoms=%s' %
          (result.get('valid'), result.get('n_atoms')))

    # 9. render is PARKED (Qt 6.11 capture-breaks-resize regression): the endpoint
    # must NOT be registered, so an agent cannot trigger it.
    if 'render' in methods:
        _fail("'render' should be parked (ENABLE_RENDER False) but is in the manifest")
    status, result = _post(port, token, 'render', width=200, height=150)
    if isinstance(result, dict) and 'error' in result and 'No registered server method' in result['error']:
        print('PASS: render correctly parked (endpoint not registered)')
    else:
        _fail('expected render to be unregistered, got: %s' % result)

    # 10. async job pattern: start a job, poll to completion, read result.
    status, result = _post(port, token, 'job_start', name='isolde status', args={})
    if status != 200 or 'job_id' not in result:
        _fail('job_start failed: %s %s' % (status, result))
    job_id = result['job_id']
    final = None
    for _ in range(100):
        status, jr = _post(port, token, 'job_poll', job_id=job_id)
        if jr.get('state') in ('done', 'failed'):
            final = jr
            break
        time.sleep(0.05)
    if final is None:
        _fail('job did not finish')
    if final['state'] != 'done':
        _fail('job failed: %s' % final)
    if 'result' not in final or 'result' not in final['result']:
        _fail('job result malformed: %s' % final)
    print('PASS: async job start/poll/result completed (state=%s)' % final['state'])

    rest_cmd.stop_server(session)
    print('ALL PASS')


try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
