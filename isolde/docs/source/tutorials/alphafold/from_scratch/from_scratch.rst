.. _alphafold_top_down_tutorial:

Top-Down Modelling with AlphaFold
=================================

**(NOTE: Most links on this page will only work correctly when the page is
loaded in ChimeraX's help viewer. You will also need to be connected to the
Internet. Please close ISOLDE and any open models before starting this
tutorial.)**

**(The instructions in the tutorial below assume you are using a wired mouse
with a scroll wheel doubling as the middle mouse button. While everything
should also work well on touchpads in Windows and Linux, support for Apple's
multi-touch touchpad is a work in progress. Known issues with the latter are
that clipping planes will not update when zooming, and recontouring of maps is
not possible.)**

Tutorial: Fitting emd_9118 using 6mhz and AlphaFold
---------------------------------------------------

.. toctree::
    :maxdepth: 2

**NOTE: This tutorial assumes that you are already familiar with the basic
operation of ISOLDE. If this is your first time using it, it is highly
recommended that you first work through the** :ref:`isolde_cryo_intro_tutorial`
**tutorial before attempting this one.**

In the :ref:`_bulk_fitting_tutorial` tutorial we used 6mhz_, the 4.1 Angstrom
structure of the ATP-bound *E. coli* LptB2FGC complex, to fit EMD-9118_, the 4.0
Angstrom map of the ATP-free, lipid-bound form. While quite achievable with
patience, this approach suffered from two key burdens: (1) unsurprisingly given
its resolution, the starting model contained numerous imperfections leading to
challenges fitting the details; and (2) like most structures, the starting model
was incomplete.

 .. _6mhz: https://www.rcsb.org/structure/6mhz
 .. _EMD-9118: https://www.emdataresource.org/EMD-9118

The advent of `AlphaFold 2`_ (AF2) (and in particular their `online database`_
of precomputed models) has fundamentally and irrevocably changed the face of
experimental model building. As demonstrated in `their paper`__ and others
(including `our own`_), in most cases the accuracy of AF2 predictions for
well-structured regions approaches that of high-resolution experimental
structures, and dramatically improves upon most existing low-resolution
structures (at least, those for which there was no high-resolution reference
model at the time of building and deposition). Even better, AF2 turns out to be
very, very good at predicting just *how* accurate its predictions are, both at a
per-residue level via its predicted local distance difference test (pLDDT)
scoring, and in its estimate of the distance error between each pair of residues
via the predicted aligned error (PAE) matrix (more on those measures later).
This promises to dramatically ease the "average" model-building task: where
previously it was commonly necessary to build residue-by-residue into density by
some combination of automatic and manual methods, in most cases it is now both
easier and less error-prone to take a "top-down" approach, starting from
near-complete models comprising the high-confidence regions of AF2 predictions
for each chain and flexibly fitting into your new map.

.. _online database: https://alphafold.ebi.ac.uk/
.. _AlphaFold 2: https://www.nature.com/articles/s41586-021-03819-2
__ `AlphaFold 2`_
.. _our own: https://www.biorxiv.org/content/10.1101/2021.09.26.461876v1

In this tutorial, we're going to revisit the 6mhz/EMD-9118 case, but this time
we're going to replace 6mhz with the AF2 predictions for each of its chains 
first. As you'll see, this makes the modelling task substantially easier - not
just because the starting model is better, but because ISOLDE's restraint schemes
can take advantage of AF2's confidence scores to more intelligently adjust their 
properties to avoid over-restraining to low-confidence predictions.

To avoid re-treading old ground, we'll start at the point of the bulk fitting 
tutorial where chains A and G are already rigid-body fitted into the map. 
`Click here`_ to set that up.

.. _Click here: `cxcmd:open\ 6mhz;open\ 9118\ from\ emdb;hide;cartoon;color\ bych;
    color\ byhet;volume\ \#2\ level\ 0.1; 
    view\ initial\ \#1;\ view\ matrix\ model\ \#1,0.198,-0.674,0.699,86.9,0.821,0.511,0.260,-63.7,-0.534,0.522,0.655,53.6;
    fitmap\ \#1/A,G\ inMap\ \#2`