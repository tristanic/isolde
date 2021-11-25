.. _remote_control_cmd:

Remote Control Options
----------------------

.. contents::
    :local:

isolde remote rest start
========================

Syntax: isolde remote rest start [**port** *integer*]

Start a background REST HTTP server to listen for commands. If no port is
specified, a random available port will be chosen and reported to the log.

isolde remote rest stop
=======================

If the REST server is running, stop it.

isolde remote rest info
=======================

Report the host address and port of the server to the log.

isolde remote xmlrpc
====================

