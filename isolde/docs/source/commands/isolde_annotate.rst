.. _isolde_annotate:

Live validation markup (isolde annotate)
========================================

The ``isolde annotate`` commands add ISOLDE's live, real-time validation
*markup* directly onto a model, updating with every change in coordinates. They
replace the former top-level ``rama`` / ``rota`` markup commands (now
deprecated). For *text* reports of outliers rather than on-model markup, use the
companion :ref:`isolde validate <isolde>` commands instead.

The bare ``isolde annotate`` command (with no subcommand) lists the available
markup types.

.. _annotate_ramachandran:

isolde annotate ramachandran
----------------------------

Syntax: isolde annotate ramachandran [*structures*] [**showFavored** *true/false* (true)]

Add live Ramachandran-validation markup to each of the given *structures* (or to
every atomic structure in the session if none are given). Each residue gets a
marker on its C-alpha atom, coloured by conformational probability (green =
favoured, shading through yellow to hot pink), and *cis*/twisted peptide bonds
are flagged with a pseudo-planar surface.

*showFavored*: if false, markers are drawn only on residues outside the favoured
regions (prior probability < 2%).

To remove the markup::

    isolde annotate ramachandran stop [*structures*]

.. _annotate_rotamers:

isolde annotate rotamers
------------------------

Syntax: isolde annotate rotamers [*structures*]

Add live rotamer-validation markup to each of the given *structures* (or to
every atomic structure in the session if none are given).

To remove the markup::

    isolde annotate rotamers stop [*structures*]
