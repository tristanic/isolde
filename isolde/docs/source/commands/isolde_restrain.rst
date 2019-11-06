.. _adaptive_distance_restraints_cmd:

Adaptive Distance Restraints
----------------------------

.. contents::
    :local:

.. _isolde_restrain_distances_cmd:

isolde restrain distances
=========================

Syntax: isolde restrain distances *atoms* [**templateAtoms** *atoms*]
[**protein** *true/false* (true)]
[**nucleic** *true/false* (true)]
[**customAtomNames** *list of names*]
[**distanceCutoff** *number* (8.0)]
[**alignmentCutoff** *number* (5.0)]
[**wellHalfWidth** *number* (0.05)]
[**kappa** *number* (5.0)]
[**tolerance** *number* (0.025)]
[**fallOff** *number* (4.0)]

Creates a "web" of adaptive distance restraints between nearby atoms,
restraining them either to their current geometry or to that of a template.
The *atoms* and *templateAtoms* arguments will be promoted to complete residues.
If *templateAtoms* is not specified, the template is the model itself.

If no template is specified, then *atoms* may simply be any selection from your
model. **(NOTE: if the selection string includes commas, it will need to be
quoted - e.g. isolde restrain distances "#1/A,B")**. If a template *is*
specified, then *atoms* and *templateAtoms* should be matched comma-separated
lists of selections, where the residues within each selection come from a
single chain, e.g. **isolde restrain distances #1/A,#1,B,#1/C:1-50
templateAtoms #2/A,#2/B,#2/C**. The individual chain selections need not have
the same number of residues, but they should of course specify closely related
sequences.

The algorithm ISOLDE uses to determine which atoms to restrain is as follows:

1. The set of atom names to use for restraints is constructed from the
   *protein*, *nucleic* and *customAtomNames* arguments. If *protein* is
   true, the following atom names are added to the list of candidates:
   - CA, CB, CG, CG1, OG, OG1

   If *nucleic* is true, these atoms are added:
   - OP1, OP2, C4', C2', O2, O4, N4, N2, O6, N1, N6

   You may, if you wish, further extend these defaults with a comma-separated
   list of other (non-hydrogen) atom names with the customAtomNames argument,
   but this should rarely be necessary.
2. A sequence alignment is performed for each pair of chains, yielding lists
   of corresponding residues.
3. All residue pairs found in the sequence alignment are merged into a single
   "super-alignment"
4. The largest rigid-body alignment of these residues' "principal atoms"
   (CA for amino acids, C4' for nucleotides) is found for which no principal
   atom is more than *alignmentCutoff* from its counterpart.
5. Within the aligned residues, matching lists are created consisting of
   atoms with allowed names present in both model and template. For each
   atom in the "template" list, all other atoms in the "template" list
   within *distanceCutoff* of the atom (excluding atoms from the same
   residue) are used to set the target distance for a restraint between
   the corresponding atoms in the "model" list.
6. Steps 4 and 5 are repeated on the residues left over from the previous
   round, until no further alignments of at least 3 residues are found.

The remaining arguments relate to the form of the restraint scheme, which
requires some explanation. The functional form is as follows:

.. math::
    E = \kappa *
    \begin{cases}
        0, & \text{if}\ enabled < 0.5 \text{ or}\ |r-r_0| < \tau \\
        1/2 (\frac{r-\rho}{c})^2, & \text{if}\ \alpha = 2 \\
        ln(\frac{1}{2} (\frac{r-\rho}{c})^2 + 1), & \text{if}\ \alpha = 0 \\
        \frac{|2-\alpha|}{\alpha} ((\frac{ (\frac{r-\rho}{c})^2 }{|2-\alpha|} + 1)^\frac{\alpha}{2} - 1), & \text{otherwise}
    \end{cases}

where

.. math::
    \rho =
    \begin{cases}
        r-\tau, & \text{if}\ (r-r_0) < -\tau \\
        r+\tau, & \text{if}\ (r-r_0) > \tau
    \end{cases}

... leading to energy potentials that look like this:

.. figure:: images/adaptive_energy_function.png

To interpret this, keep in mind that the force applied to a given atom is
proportional to the *derivative* (that is, the slope) of the energy with
respect to distance. In effect, for a pair of atoms that is close to the target
distance the restraint will act as a simple harmonic spring (with a
"flat-bottom" range defined by *tolerance* over which no force is applied).
Once the interatomic distance deviates from the target by more than
(*wellHalfWidth* + *tolerance*) * target_distance the energy profile begins to
flatten out at a rate specified by :math:`\alpha`.

The value of :math:`\alpha` for a given restraint is determined by the
combination of *fallOff* and the target distance, based on the idea that larger
distances are inherently less certain. Specifically,
:math:`\alpha = -2 -\text{fallOff} ln(\text{target})`. Small (or negative
- **not** generally recommended) values of *fallOff* will cause the restraints
to stay quite strong with increasing distances (like a ball rolling into a
funnel), whereas large values cause them to quickly flatten off (like a golf
ball rolling into its cup - only falling into place when *very* close to the
target).

Finally, the parameter *kappa* sets the overall strength of the restraints.
The effective spring constant for a given restraint within its "harmonic" range
is :math:`k=\frac{\kappa}{(\text{wellHalfWidth}*\text{(target distance)})^2}`
:math:`kJ mol^{-1} nm^{-2}`.

isolde restrain single distance
===============================

Syntax: isolde restrain single distance *atoms* *minDist* *maxDist*
[**strength** *number* (20)]
[**wellHalfWidth** *number* ((*minDist*+*maxDist*)/10)]
[**confidence** *number* (-2)]

Restrain the distance between a single pair of atoms. The arguments *minDist*
and *maxDist* specify the range of distances between which no force will be
applied. In other words, the target distance becomes (*minDist*+*maxDist*)/2,
and the tolerance becomes (*maxDist*-*minDist*)/2.

* *strength*: Corresponds directly to :math:`\kappa` in the energy formula.
* *wellHalfWidth*: Corresponds directly to *c* (the half-width of the harmonic
  region) in the energy formula. If not specifies, it defaults to 1/5 of the
  target distance.
* *confidence*: Corresponds to :math:`\alpha` in the energy formula.

isolde release distances
========================

Syntax: isolde release distances *atoms* [**internalOnly** *true/false* (false)]
[**externalOnly** *true/false* (false)] [**longerThan** *number*]
[**strainedOnly** *true/false* (false)] [**stretchLimit** *number* (1.2)]
[**compressionLimit** *number* (0.8)]

Release a selection of adaptive distance restraints. *(NOTE: released restraints
cannot currently be reinstated, but may be re-created using the "isolde restrain
distances" command)*

Calling *isolde restrain distances <selection>* with no other arguments will
simply release all restraints involving any of the specified atoms (including
restraints to atoms outside the selection). The remaining arguments allow fine-
tuning of the selection to release:

* *internalOnly* (**incompatible with externalOnly**): if true, only those
  restraints for which both atoms are within the selection will be
  released.
* *externalOnly* (**incompatible with internalOnly**): if true, only those
  restraints connecting atoms within the selection to those outside will be
  released.
* *longerThan*: a value in Angstroms. If specified, only restraints with
  target distances larger than this value will be released.
* *strainedOnly*: if true, only restraints with (length/target) larger than
  *stretchLimit* or smaller than *compressionLimit* will be released.
* *stretchLimit* (**ignored unless strainedOnly is true**): ratio of current
  distance to target distance above which restraints will be released.
* *compressionLimit* (**ignored unless strainedOnly is true**): ratio of
  current distance to target distance below which restraints will be
  released.

isolde adjust distances
=======================

Syntax: isolde adjust distances *atoms* [**internalOnly** *true/false* (false)]
[**externalOnly** *true/false* (false)] [**kappa** *number*]
[**wellHalfWidth** *number*] [**tolerance** *number*] [**fallOff** *number*]
[**displayThreshold** *number*]

Adjust the strength and/or display properties of a set of adaptive distance
restraints.

* *internalOnly* (**incompatible with externalOnly**): if true, only those
  restraints for which both atoms are within the selection will be
  released.
* *externalOnly* (**incompatible with internalOnly**): if true, only those
  restraints connecting atoms within the selection to those outside will be
  released.
* *kappa*: see :ref:`isolde_restrain_distances_cmd`
* *wellHalfWidth*: see :ref:`isolde_restrain_distances_cmd`
* *tolerance*: see :ref:`isolde_restrain_distances_cmd`
* *fallOff*: see :ref:`isolde_restrain_distances_cmd`
* *displayThreshold*: deviation from target distance (expressed as a fraction
  of the target distance) below which a given restraint will not be shown. For
  example, to only show restraints deviating more than 10% from their targets,
  set *displayThreshold* to 0.1. To show all restraints, set displayThreshold to
  0.
