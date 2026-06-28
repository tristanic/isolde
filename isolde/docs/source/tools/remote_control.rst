.. _remote_agent_control:

Remote control and AI agents (MCP)
==================================

.. toctree::
    :maxdepth: 2

.. contents::
    :local:
    :depth: 2

What this is
------------

ISOLDE can be driven by an external program — including a large language model
"agent" such as Claude — through the `Model Context Protocol
<https://modelcontextprotocol.io>`_ (MCP). The agent calls ISOLDE operations as
typed, discoverable tools: it can look at your model, run validation, start and
steer simulations, impose restraints, refine B-factors, and so on.

Crucially, it does this **alongside you, in your live ChimeraX/ISOLDE GUI** — not
in a hidden headless session. You see every map, every simulation frame, and
(see :ref:`agent_transparency`) a log line for every change the agent makes. The
GUI is the safety net: you remain in control and can pause, undo, or veto at any
time.

Under the hood the agent talks to ISOLDE over a small authenticated REST server
bound to ``localhost`` (command reference: :ref:`remote_control_cmd`). A separate
**MCP server** process translates between the agent's MCP calls and that REST
API.

.. _agent_philosophy_summary:

The philosophy behind it
------------------------

ISOLDE is opinionated about *how* models should be built, and the agent
interface bakes that in so an automated client does not quietly do the wrong
thing. In short:

- **The map is the evidence; the model is a hypothesis.** The goal is a model
  that explains the density with good physics — not one tuned to make a score
  look good.
- **"Non-outlier" is not the same as "correct".** A clean clashscore or
  Ramachandran plot does not prove a region is right; good scores can hide real
  errors.
- **Small errors are the dangerous ones.** Sub-Ångström switches, *cis*-peptides
  and flipped ligands pass static validation but can change the science.
- **Never force-fit.** Cranking up map weighting to chase fit wrecks geometry;
  the aim is a low-energy model that *also* fits.
- **Trust the physics; let it settle.** Express intent (nudge an atom toward
  density) and let the simulation relax around it, rather than hard-setting
  coordinates.
- **You are the gold standard.** "Human eyes should see each residue in density
  at least once." The agent works with you, not instead of you.

The agent receives these principles automatically when it connects (as MCP
*instructions* and an ``isolde://philosophy`` resource), and the live runtime
copy of that text is the single source of truth. The summary here is for human
readers.

Activating the REST server
--------------------------

In the ChimeraX command line::

    isolde remote rest start

ISOLDE picks a free port, mints an authentication token, and prints both to the
log. To see them again at any time::

    isolde remote rest info

You will need the **port** and **token** to connect a client. Stop the server
with ``isolde remote rest stop``. See :ref:`remote_control_cmd` for the full
command reference.

.. warning::

    Do **not** pass ``allowRun true`` unless you know you need it: it enables a
    method that runs arbitrary ChimeraX commands. With the default the agent is
    limited to ISOLDE's typed, opt-in tool set.

Running the MCP server
----------------------

The MCP server is a small standalone process that ships inside the ISOLDE
bundle, at ``chimerax/isolde/mcp/server.py`` under your ChimeraX installation. It
depends only on the ``mcp`` Python package and the standard library, so it can run
in any Python 3 environment::

    pip install mcp

It reads the REST connection details from environment variables and speaks MCP
over stdio:

.. code-block:: bash

    export ISOLDE_REST_HOST=localhost          # default
    export ISOLDE_REST_PORT=<port from 'isolde remote rest info'>
    export ISOLDE_REST_TOKEN=<token from 'isolde remote rest info'>
    python /path/to/chimerax/isolde/mcp/server.py

(On Windows PowerShell, use ``$env:ISOLDE_REST_PORT="..."`` etc.)

Connecting a client
--------------------

Most MCP clients (for example Claude Desktop or Claude Code) are configured by
pointing them at the command that launches the server, with the connection
details supplied as environment variables. A typical stdio entry looks like:

.. code-block:: json

    {
      "mcpServers": {
        "isolde": {
          "command": "python",
          "args": ["/path/to/chimerax/isolde/mcp/server.py"],
          "env": {
            "ISOLDE_REST_HOST": "localhost",
            "ISOLDE_REST_PORT": "<port>",
            "ISOLDE_REST_TOKEN": "<token>"
          }
        }
      }
    }

.. note::

    The port and token are **regenerated every time you run** ``isolde remote
    rest start``. After restarting the REST server, refresh the ``ISOLDE_REST_PORT``
    and ``ISOLDE_REST_TOKEN`` values (re-read them with ``isolde remote rest
    info``) and restart the MCP server / reconnect the client. Consult your MCP
    client's own documentation for exactly where its server configuration lives.

What the agent can see and do
-----------------------------

On connecting, the MCP server advertises:

- **Perception / query helpers** — read-only structured views an eyeless agent
  needs: a model inventory, a model description (chains, sequences, ligands,
  chain breaks), atom-spec resolution, clustered "problem zones" (where to look
  first), restraint and per-residue detail, density-fit metrics (Rwork/Rfree and
  MDFF state), and B-factor outliers.
- **Command tools** — one tool per ISOLDE command that has been opted in as
  ``agent_safe``: simulation control, manipulations (peptide flips, etc.),
  validation, restraints, navigation and B-factor/occupancy refinement. Exposing
  a new command is a one-line change, so this set grows over time.
- **Action helpers** — ``tug`` (a gentle, transient spring that nudges atoms
  toward a target and lets the simulation relax around them) and conservative
  control of MDFF map coupling.

Two safeguards run on top of these:

- **Witnesses.** Mutating commands return a *witness* — a verifiable statement of
  what actually changed (atoms moved, restraints added, the simulation stepping).
  This lets the agent (and you) tell real work from a command that "succeeded"
  without doing anything.
- **Soft guardrails.** When a call looks like it violates a principle above (for
  example raising MDFF coupling to chase fit, or an unusually large single move),
  the result carries an advisory ``philosophy_flag``. These *flag*, they never
  block — the human and the GUI remain the safety net.

How an agent should approach a model
------------------------------------

ISOLDE's :doc:`tutorials <../tutorials/isolde>` converge on a consistent way of
working, and an agent should follow the same loop. This is the human-readable
mirror of the workflow the agent receives at runtime in the
``isolde://philosophy`` resource (the canonical copy):

#. **Prepare for simulation** — add hydrogens, remove alternate conformations, and
   resolve unparametrised residues, so every residue is known to the force field.
#. **Settle once** — run a whole-model simulation to relieve clashes and let obvious
   errors fix themselves. You typically simulate the whole model only twice: once at
   the start to relieve clashes, once at the end to settle.
#. **Triage** — use the problem-zones query to find the worst clusters and fix the
   biggest first; many downstream issues cascade away once the root cause is fixed.
#. **Fix locally** — peptide flips, rotamer corrections, register shifts; tug in
   *short bursts* to help atoms into density, then let the dynamics settle.
#. **Reconcile restraints with the map** — where reference/AlphaFold restraints fight
   strong density (a cluster of strained restraints plus poor fit), the reference is
   probably wrong: release those restraints and let the density win.
#. **Final settle** — a whole-model simulation with the temperature dropped toward
   zero, to bring local geometry to equilibrium before writing coordinates.
#. **Inspect everything** — step through every residue in context with its density;
   human eyes remain the gold standard.
#. **Refine externally** — ``isolde brefine``/``brsr`` handle B-factors and
   occupancies; a final fine-geometry polish in a dedicated package may still help.

Rules of thumb
~~~~~~~~~~~~~~~

- Never add or remove atoms while a simulation is running.
- Don't chase fit by raising MDFF coupling; watch for over-fitting symptoms (backbone
  twisting out of stable secondary structure, sidechains pulled too aggressively into
  the map).
- Tug in short bursts to *help* the model into place; one hard pull force-fits and
  wrecks local geometry.
- Severe clashes (atoms passing through each other) won't resolve by naive
  minimisation — separate them gently (soft-core nonbonded potentials plus local
  restraints over a region big enough to let the pieces slide apart) rather than
  cranking up forces.
- Validation is a hint, not a verdict: not every outlier is wrong, and not being an
  outlier doesn't make something correct. A peptide twist beyond ~30° is effectively
  never real; a non-proline *cis* bond is rare (~3 in 10,000) and, when genuine,
  well-resolved and functionally interesting; ~5% of prolines are legitimately *cis*.
- AlphaFold / reference models: residues with pLDDT < 50 are essentially junk (don't
  interpret those coordinates), and only residue pairs with low predicted error
  (PAE ≲ 4 Å) are worth restraining — a confident prediction can still be wrong where
  the map disagrees.
- For out-of-register stretches, use the register-shifter rather than a tug-of-war
  between distance restraints.

.. _agent_transparency:

Seeing what the agent does
--------------------------

Because you share the live session, every **mutating** agent action is echoed to
the ChimeraX log as a concise ``[agent] ...`` line, followed by its witnessed
outcome — for example::

    [agent] isolde restrain torsions #1.2  → +815 restraints

Read-only queries (validation, status, model inspection) stay silent to avoid
flooding the log. High-severity guardrail flags are also surfaced in the log, so
a questionable move (such as a coupling increase) is hard to miss.

Limits and safety
-----------------

- The REST server binds to ``localhost`` only and requires the bearer token on
  every request.
- The free-text ``run`` method is **off by default** (``allowRun``); without it
  the agent is confined to the typed, opt-in tool set.
- **Agent vision (rendering the view back to the agent as an image) is currently
  disabled**, pending a Qt regression that affects offscreen capture. Until it is
  re-enabled, the agent works from structured queries and your shared view rather
  than its own snapshots.
- Not all chemistry is covered: real-time validation focuses on protein
  torsion-based metrics, and the force field covers proteins, nucleic acids,
  glycans and a large set of common ligands but not arbitrary novel chemistry.
  Outside that coverage the agent should flag the gap rather than guess.

See also
--------

- :ref:`remote_control_cmd` — the ``isolde remote rest`` command reference.
- :ref:`isolde_gui` — the interactive GUI the agent shares with you.
