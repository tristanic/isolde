Live rotamer validation
=======================

Each protein sidechain with at least one rotatable bond (that is, excluding
alanine and glycine) has a fairly limited range of stable (that is,
energetically-favourable) conformations known as *rotamers*. Statistical
analyses of very large numbers of well-resolved residues from high-resolution
crystal structures has left us with detailed joint probability distributions
for the torsion angles around each rotatable sidechain bond (known as *chi*
angles). These distributions are contained in the `MolProbity ultimate rotamer library`_
and are used by ISOLDE to provide real-time rotamer validation as a model's
coordinates change. Use of this validation is not restricted to ISOLDE, though:
the *rota* command can be used to add live rotamer validation markup to any
atomic model in ChimeraX.

.. _MolProbity ultimate rotamer library: https://onlinelibrary.wiley.com/doi/full/10.1002/prot.25039

The live validation markup looks like this:

.. figure:: images/rota_markup.jpg

Each residue whose conformation has a prior probability of less than 2% will
have the pictured exclamation mark/spiral motif parallel to its CA-CB bond. Both
scale and colour of the indicator change with the severity of the outlier: for a
residue with a probability close to 2% the indicator will be small and yellow;
at a probability of 0.05% the indicator reaches its maximum size and strongest
red/pink colour. The indicators are automatically updated whenever atoms move.

Syntax: rota [*structures*] [**report** *true/false* (false)]


*structures*: if provided, a validator will be added to each of the specified
structures. Otherwise, one will be added to each structure in the ChimeraX
session.

*report*: if true, a report listing all residues with probabilities less than
0.02 will be printed to the log.

.. _`stop`:

rota stop
---------

Syntax: rota stop [*structures*]
Alias: ~rota [*structures*]

Remove any live rotamer validators from the specified structures. If *structures*
is not specified, all rotamer validators will be removed from the current
ChimeraX session. This command is precisely equivalent to simply
closing the "Rotamer Validation" model in the Models panel.
