.. _remote_control_cmd:

Remote Control Options
----------------------

.. contents::
    :local:

ISOLDE can expose a running session to other local programs through a small
**authenticated REST server** bound to ``localhost``. This is the substrate the
:ref:`MCP / AI-agent interface <remote_agent_control>` rides on, but it is a
general-purpose HTTP API you can also drive from your own scripts. For the
friendly, "what is this and how do I use it" overview, see
:ref:`remote_agent_control`; this page is the command reference.

Every request must carry the bearer token minted when the server starts (see
``isolde remote rest start`` / ``info`` below); requests without it are rejected
with HTTP 401.

isolde remote rest start
========================

Syntax: isolde remote rest start [**port** *integer*] [**allowRun** *true/false* (false)]

Start a background REST HTTP server, listening on ``localhost`` only. If no
**port** is given, a random free port is chosen. On start the server mints a
random **authentication token** and prints it to the log; every request must
supply it as an ``Authorization: Bearer <token>`` header. The port and token are
also available later via ``isolde remote rest info``.

**allowRun** (default ``false``) enables the free-text ``run`` method, which
executes arbitrary ChimeraX commands sent by a client.

.. warning::

    ``allowRun true`` lets any local program that has the token run *arbitrary
    ChimeraX commands* in your session. Leave it off unless you specifically need
    it, and only enable it for trusted local automation. With the default
    (``false``) the agent surface is restricted to the typed, opt-in set of
    ``agent_safe`` commands and helpers.

isolde remote rest info
=======================

Report, to the log, whether the server is running and — if so — the host, the
port, the **authentication token**, and whether the free-text ``run`` method is
enabled. A client (for example the MCP server) reads the port and token from
here to connect.

isolde remote rest stop
=======================

If the REST server is running, stop it and release its port. (The token from the
stopped server is no longer valid; starting again mints a fresh one.)

isolde remote xmlrpc
====================
