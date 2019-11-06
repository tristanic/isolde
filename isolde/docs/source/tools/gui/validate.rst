Monitoring the quality of your model
====================================

.. toctree::
    :maxdepth: 2

.. contents::
    :local:

While ISOLDE's live annotations help you see at a glance where the problems are
in whatever part of the model you happen to be looking at, it is of course still
important to have access to a quick summary of *all* its known (potential)
problems. This is handled on the `Validate` tab on the ISOLDE panel. If you
click there, you'll intially be greeted by this rather boring view:

.. figure:: images/validate_tab_all_closed.png
    :alt: Initial view of the Validate tab

    Expect the number of options to grow as ISOLDE develops further

Let's go through and look at what each of these buttons does.

The Ramachandran plot
---------------------

Clicking the **Show** button on the *Ramachandran Plot* bar will give you a view
that looks something like this:

.. figure:: images/validate_tab_rama_plot.png
    :alt: The Ramachandran plot

    The Ramachandran plot for "general" amino acids

In brief, the Ramachandran plot plots the angle of the *phi* (φ) dihedral
"entering" the residue against that of the *psi* (ψ) dihedral "exiting" it.

.. figure:: images/phi_dihedral.png
    :alt: The phi dihedral

    The four atoms making up the *phi* (φ) dihedral.

.. figure:: images/psi_dihedral.png
    :alt: The psi dihedral

    The four atoms making up the *psi* (ψ) dihedral.

Since atoms are physical objects that *really* don't like to overlap with each
other, only some combinations of φ and ψ are possible (and only a much smaller
proportion of the range of possible combinations are *favourable*). That is
what the Ramachandran plot is telling you.

The yellow contour line represents  the 2% cut-off: 98% of all well-resolved
residues in very high-resolution  structures are found within these contours.
Unsurprisingly, they mostly  represent the well-defined secondary structures:
the big yellow blob in the top left corresonds to α-helices, while the slightly
smaller middle-left blob  corresponds to β-strands. Residues within these
contours are said to be in  "favoured" conformations.

The purple line, meanwhile, represents the boundary of the far more generous
"allowed" region. The precise value of this contour varies slightly depending on
the Ramachandran category, but for the general case only about 5 in 10,000
residues should be found outside it, in what is known as "outlier" space.

Note that as well as the given cut-off contours, the background is shaded
according to probability as well. A residue outside the contour but still on a
grey zone may still be correct (but should be very carefully checked), but one
sitting on a white background is almost certainly wrong.

As you can see from the example image above, the size and colour of each point
scales with its probability, making outliers very difficult to ignore. Clicking
on any point in the plot will take you to the corresponding residue in the main
window and select it, so (if you don't already have a simulation running),
clicking the play button will immediately start a localised simulation allowing
you to work on it.

The drop-down menu below the plot allows you to limit its scope to only the
residues currently selected in the main *ChimeraX* window. This can be useful
if, for example, you wish to consider only a single chain. While a simulation is
running, the plot will be limited to only the mobile residues. While visible,
the Ramachandran plot updates automatically whenever the atomic coordinates
change.

Peptide Bond Geometry validation
--------------------------------

.. figure:: images/validate_tab_peptide_geometry.png
    :alt: Peptide bond geometry validation

    Tabulated list of questionable peptide bonds

    Red: non-proline *cis* peptide bond
    Green: proline *cis* peptide bond
    Yellow: twisted peptide bond (more than 30° from planar)

Clicking on any entry in the list will select the residue in question and take
you to it in the main view.

In a peptide bond, the valence electrons are partially delocalised (think of it
as the carbon-nitrogen bond spending some fraction of its time as a double
bond), giving it a strong bias towards planarity. While in some rare, tightly
constrained conditions real peptide bonds have been observed to be twisted up to
30°, any twist beyond that is effectively unheard of.

The question of *trans* vs *cis* peptide bonds is less obvious. The vast
majority of all peptide bonds are found in *trans* (with the amide hydrogen
pointing in the opposite direction to the carbonyl oxygen):

.. image:: images/transpep.png
    :alt: a normal trans peptide

This occurs due to steric considerations: swapping to *cis* replaces the amide
hydrogen with the next alpha carbon, and carbon atoms are simply too big to
pack neatly together in this way. Forcing them to do so puts strain on the
intervening bonds, which requires energy. For this reason, for all residues
other than proline *cis* peptide bonds are vanishingly rare at about 3 per
10,000 residues - and when they *do* occur are strongly stabilised by
interactions with their surroundings, usually well-resolved, and almost always
functionally interesting.

.. figure:: images/cispep_loop.png
    :alt: a series of cis peptide bonds

    This should never happen.

    A loop like this, with four non-proline *cis* peptide bonds in a row and no
    stabilising influences whatsoever, is for all intents and purposes
    impossible.

.. figure:: images/real_cis_bond.png
    :alt: a real cis peptide bond

    This, on the other hand, is real.

    This non-proline *cis* bond is found in tissue transglutaminase (*TGM2*, see
    PDB ID 2q3z), and is stabilised by a disulphide bond between directly
    adjacent cysteine residues. It is part of a redox-switching mechanism: under
    certain circumstances the disulphide is reduced by thioredoxin, allowing
    the peptide bond to flip to *trans* and thereby shifting an inhibitory loop
    to activate the enzyme. It is also clearly resolved in the electron density.

If your site in question does not meet the above criteria, then in general
you should assume it to be *trans* (unless, of course, you have outside
evidence such as a higher-resolution structure of a closely-related protein).

Proline is a special case. Unlike the other amino acids, proline does not have
an amide hydrogen - instead the amide nitrogen links to the backbone delta
carbon. With both *cis* and *trans* conformers putting a carbon atom adjacent to
the alpha carbon, the drive towards the *trans* conformation is nowhere near as
strong, and about 5% of all proline residues appear in *cis*. With this in mind,
*cis* prolines are coloured in green to reflect their higher likelihood (but
keep in mind that they do still need to agree with the density map, as in the
below case).

.. figure:: images/cis_pro.png
    :alt: a cis-proline bond

    A happy-looking *cis*-proline.

Rotamer Validation
------------------
.. figure:: images/validate_tab_rotamers.png
    :alt: the rotamer validation table

All non-favoured rotamers are listed in the table, ordered from worst to best.
As for the other validation tools, clicking on an entry will select and take you
to the offending residue.

Clashes
-------

The behaviour of this widget is subtly different depending on whether a
simulation is currently running. Outside of simulations, clashes are determined
based on simple pairwise distances:

.. figure:: images/clash_table_no_sim.png
    :alt: clash widget with no simulation running

While a simulation is running, the table is instead populated by a list of atoms
experiencing large net forces:

.. figure:: images/clash_table_during_sim.png
    :alt: clash widget during a simulation

In either case, clicking on an entry in the table will focus the view on the
offending atom(s). In general, during a simulation the clash table will only be
populated (a) if updated before minimisation is complete; or (b) if the
minimiser is unable to resolve one or more severe clashes (see
:ref:`dealing-with-clashes`).
