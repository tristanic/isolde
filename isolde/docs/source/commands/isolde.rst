The ISOLDE Command Line
-----------------------

.. toctree::
    :maxdepth: 2

    isolde_restrain
    isolde_remote

While the primary mode of control of ISOLDE is via :ref:`isolde_gui`, if you
prefer you can also perform many tasks via the command line, such as launching
the GUI (:ref:`start`), starting, stopping and pausing simulations (:ref:`sim`),
and performing basic manipulations (:ref:`pepflip` and :ref:`cisflip`). Some
functions are currently available *only* through the command line - in
particular, :ref:`adaptive_restraint_schemes`.


.. _tutorial:

isolde tutorial
===============

Brings up the ISOLDE :ref:`isolde_tutorials` help page, providing interactive
case studies for working on models in ISOLDE.

.. _demo:

isolde demo
===========

Syntax: isolde demo *demoName* [**modelOnly** *true/false*] [**startIsolde** *true/false*]

Primarily intended for use with the ISOLDE tutorials. Loads saved atomic
coordinates for use in rebuilding.

*demoName*: either "cryo_em" or "crystal"

*modelOnly*: if true, load only the model and not the electron density map

*startIsolde*: if true, the ISOLDE GUI will be started.


.. _start:

isolde start
============

Brings up :ref:`isolde_gui`. The first model in the list of open models will be
automatically prepared with ISOLDE's default visualisation mode. This also sets
the display camera to orthographic projection and lighting to simple. It is
strongly recommended that you do not change these settings during your ISOLDE
session.

.. _set:

isolde set
==========

Syntax: isolde set [**timeStepsPerGuiUpdate** *integer*]
[**temperature** *number*] [**gpuDeviceIndex** *integer*]

Change various ISOLDE global settings.

isolde reset forcefield
=======================

Delete cached forcefield files to reset to the as-installed state.

.. _report:

isolde report
=============

Sytax: isolde report [**true|false** (true)]
[**interval** *integer* (20)]

Start/stop reporting information on simulation performance (time per coordinate
update and timesteps per second) to the status bar. The optional
*interval* argument sets the number of coordinate updates to average over
before reporting. Only valid while a simulation is running, and automatically
terminates once that simulation stops.

.. _sim:

isolde sim
==========

Syntax: isolde sim *cmd* [*atoms*] [**discardTo** *discardTo*]

Start, stop or pause an interactive simulation.

isolde sim start [*atoms*]
~~~~~~~~~~~~~~~~~~~~~~~~~~

Launches the ISOLDE GUI if not already running, and starts a simulation. If no
atoms are specified, the simulation will encompass the entirety of the current
selected model (or the first open model if ISOLDE is not already running). If
*atoms* is specified, the selected atoms must come from a single model
(if this is not the current selected model, ISOLDE will automatically switch).

isolde sim stop [**discardTo** *discardTo*]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*discardTo*: one of "checkpoint" or "start"

Stop the currently running simulation, and optionally discard the results. If
*discardTo* is not specified, the current state of the model will be kept.
Otherwise, the atomic positions and states of all restraints will be reverted
either to the starting state or the last saved checkpoint.

isolde sim pause
~~~~~~~~~~~~~~~~

Pauses/resumes the current simulation.

.. _pepflip:

isolde pepflip
==============

Syntax: isolde pepflip *atoms*

Attempt to flip the peptide bond N-terminal to each selected residue, starting a
suitable simulation if none is currently running. Requires ISOLDE to already be
initialised. Residues outside the model currently selected for ISOLDE (or
outside the mobile selection if a simulation is running) will not be flipped.

.. _cisflip:

isolde cisflip
==============

Syntax: isolde cisflip *atoms*

Attempt to flip the peptide bond N-terminal to each selected residue from *cis*
to *trans* or vice versa, starting a suitable simulation if none is currently
running. Requires ISOLDE to already be initialised. Residues outside the model
currently selected for ISOLDE (or outside the mobile selection if a simulation
is running) will not be flipped.

.. _ignore:

isolde ignore
=============

Syntax: isolde ignore *residues*

Tell ISOLDE to ignore a selection of residues in future simulations. This will
not take effect until the next simulation is started. Ignored residues will
still be used for structure factor calculations, but do not take part in
simulations in any way. Atoms directly bonded to an ignored residue will be
fixed in space. This command is most useful in dealing with severe clashes that
are otherwise intractable (e.g. docked homology models with intertwined loops):

* select one of the clashing chains
* *isolde ignore sel* to ignore it
* start a simulation and remodel its counterpart into density
* stop the simulation, then *isolde ~ignore* to clear the list of ignored
  residues
* if necessary, ignore the selection you just remodelled and run a
  simulation to fix the other side of the problem area
* stop any running simulation, *isolde ~ignore*, then run a final simulation
  with everything included to resolve any remaining clashes.

.. _`~ignore`:

isolde ~ignore
==============

Syntax: isolde ~ignore *residues*

Tell ISOLDE to stop ignoring a selection of residues for simulation purposes.
This will not take effect until the next new simulation is started.

.. _`stepto`:

isolde stepto
=============

Syntax: isolde stepto [*residue or {next|prev}*] [**viewDistance** *number*]
[**interpolateFrames** *integer*] [**polymericOnly** {*TRUE|false*}]

Focus the camera on the specified residue, or the next/previous residue in the
chain if no residue is specified. *isolde stepto next* will move the camera to
the next residue in the model, while *isolde stepto prev* will move back to the
previous one. The stepper will remember the last specified direction, so
repeated calls of "isolde stepto" without arguments will keep moving along the
chain in the same direction. Providing the *viewDistance* argument will cause
the camera to zoom in/out to the specified distance; this distance will be
maintained for future calls. If the current camera position is close, the view
will slide smoothly to the new position over *interpolateFrames* frames, otherwise
it will jump directly there. If *polymericOnly* is true (default), any residues
not part of a protein or nucleic acid chain will be skipped. As for the other
arguments, this will be remembered for all future calls.

Each loaded model is given its own independent residue stepper - the settings
you make for one will not be carried over to others.

.. _`jumpto`:

isolde jumpto
=============

Syntax: isolde jumpto [*next|prev*]

Jump the residue stepper to the first residue of the next chain, or last residue
of the previous chain.

.. _`add ligand`:

isolde add ligand
=================

Syntax: isolde add ligand *{residue ID}* [*model*]
[**position** *list of three floats*] [**bfactor** *number*] [**chain** *string*]
[**distanceCutoff** *number*] [**simSettle** {*true|FALSE*}]
[**useMdTemplate** {*TRUE|false*}] [**mdTemplateName** *string*]

*NOTE: when placing ligands using this command, ISOLDE does not currently make
any attempt at a preliminary fit prior to starting simulations - it will simply
place the "ideal" coordinates specified in the template file. When adding large,
flexible ligands this will almost always lead to severe clashes with the
surroundings. In such cases, it is advisable to use* **isolde ignore ~sel** *to
exclude everything but the ligand from simulations, perform an initial fit to
the density using tugging and/or position restraints on key atoms, then use
* **isolde ~ignore** *to reinstate the rest of the model for simulations and
continue on.*

Add a ligand based on its template in the Chemical Components Dictionary. The
residue ID must match the 3-letter code for the ligand (if you don't know this,
you can search for it at http://ligand-expo.rcsb.org/). If ISOLDE is started,
then the *model* argument is optional (if not provided, the ligand will be added
to ISOLDE's currently selected model). Otherwise, you will need to explicitly
specify the model to add to. By default, the ligand will be added at the current
centre of rotation; you may, however, specify an alternative location by
providing coordinates as a comma-separated list after the *position* keyword.

By default, the b-factor and chain assigned to the residue will be determined
based on the closest atoms within *distanceCutoff* (default: 8 Angstroms) to the
site of addition, but you may explicitly specify these if you wish. If
*simSettle* is true, a local simulation will automatically be started - this is
only advisable for small, rigid molecules for which severe clashes with the
surroundings are unlikely.

If *useMdTemplate* is true, the added residue will be checked against the
corresponding molecular dynamics template (if present), and atoms will be added
or removed as needed to match (templates provided by the CCD are often not in
the protonation states most common under biological conditions). You should not
usually need to use the *mdTemplateName* argument: if you have loaded a custom
template, it will be found and used as long as its name matches the residue
name.

.. _`add water`:

isolde add water
================

Syntax: isolde add water [*model*] [**position** *list of three floats*]
[**bfactor** *number*] [**chain** *string*]
[**distanceCutoff** *number*] [**simSettle** {*TRUE|false*}]

Essentially a special case of *isolde add ligand*. The primary difference is
that *simSettle* defaults to true (that is, adding a water will automatically
start a local simulation to settle it). In addition, the default value for
*distanceCutoff* is reduced from 8.0 to 3.0 Angstroms, on the basis that it is
rarely a good idea to add a water outside of hydrogen bonding distance from the
nearest existing atom(s).

.. _`replace ligand`:

isolde replace ligand
=====================

Syntax: isolde replace ligand *residue* *newResidueName*

*(EXPERIMENTAL)*

Replace one ligand with a related one, keeping as many atoms common to both as
possible. Matching of common atoms is performed by graph matching based on
bonding between elements. Use with caution: the current implementation is not
aware of bond order nor of chirality, so attempting to replace (for example) a
D-sugar with its L-enantiomer will simply rename the residue while retaining the
D coordinates. This will be improved upon in a future release.



.. _`restrain distances`:

isolde restrain distances
=========================

See :ref:`adaptive_restraint_schemes`

.. _`release distances`:

isolde release distances
========================

See :ref:`adaptive_restraint_schemes`

.. _`adjust distances`:

isolde adjust distances
=======================

See :ref:`adaptive_restraint_schemes`


.. _`restrain torsions`:

isolde restrain torsions
========================

See :ref:`adaptive_dihedral_restraint_cmd`

.. _`adjust torsions`:

isolde adjust torsions
======================

See :ref:`adaptive_dihedral_restraint_cmd`

.. _`release torsions`:

isolde release torsions
=======================

See :ref:`adaptive_dihedral_restraint_cmd`
