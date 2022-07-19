.. _alphafold_mr_tutorial:

Starting From a Molecular Replacement Solution for a Flexible Protein
=====================================================================

**NOTE: Most links on this page will only work correctly when the page is
loaded in ChimeraX's help viewer. Unless the necessary files have been
pre-fetched you will also need to be connected to the Internet. Please close any
open models before starting this tutorial.**

**The instructions in the tutorial below assume you are using a wired mouse
with a scroll wheel doubling as the middle mouse button. While everything
should also work well on touchpads in Windows and Linux, support for Apple's
multi-touch touchpad is a work in progress. Known issues with the latter are
that clipping planes will not update when zooming, and recontouring of maps is
not possible.**

If you wish to leave this tutorial part-way through and return to it later, you
can save the state (all open models, maps, restraints etc.) to a ChimeraX session 
file via the File/Save menu option or with the command:

*save tutorial_session.cxs*

You can of course also save the working model to PDB or mmCIF - again, via the 
menu or with your choice of the below commands:

*save model.pdb #1*

*save model.cif #1*

To save without hydrogens, instead do:

*sel ~H; save {your choice of filename and extension} selectedOnly true*

The above commands will save to your current working directory. They are not 
provided as links here because that would save them to the directory containing
this tutorial file, which is probably *not* what you want.

Tutorial: Flexibly fitting an AlphaFold model to match a molecular replacement solution
---------------------------------------------------------------------------------------

*Level: Moderate*

.. toctree::
    :maxdepth: 2

**This tutorial will write a number of files to your current working directory.
To avoid clutter and confusion it is best if you create and change to a fresh
directory first.**

`Click here to choose a directory and populate it with tutorial data`__

__ cxcmd:cd\ browse;isolde\ tut\ prep\ alphafold_mr

Solving a crystallographic dataset by molecular replacement differs from fitting into
cryo-EM density in two important ways. First, there is the fact that the quality of 
a crystallographic map is dependent on both the quality and completeness of the atomic
model: too many missing ordered atoms, or two many atoms far out of position, and the 
map degrades to the point of uninterpretability. Second, the crystallographic model is 
one repeat unit in a continuous array of contacting copies - out-of-position atoms will 
commonly end up overlapping with neighbouring repeat units.

The problem is that the vast majority of proteins like to *move* at various scales, from
sidechains and loops to large-scale interdomain rearrangements. In these cases, successful
molecular replacement relies on first breaking down your search model into rigid fragments,
including trimming away mobile loops. It is common for 10-50% of protein residues to be 
discarded in this process, which can be a real shame when your search model is actually 
good quality, just in a different conformation. Historically, the standard approach has 
been to extend from the MR solution by iteratively tracing as much of the residual density
as possible (manually or via various automatic tools) interspersed with rebuilding and 
refining. Depending on the quality of the dataset this can be quite slow, painful and 
error-prone.

In this tutorial we will explore a way to short-circuit this task - instead of building out
from the MR solution, we will use it (and to start with, the initial map generated from it)
as a guide to refit the *complete* original model (in this case an AlphaFold model, but 
essentially the same procedure could be used with an existing experimental model if you 
have one). As our example, we'll be recapitulating `3now`_, the 3 Angstrom structure of 
*Drosophila* UNC-45. This largely alpha-helical structure forms a lopsided "V", and 
AlphaFold models disagree with the experimental structure primarily in the tightness of 
the central angle allowing rigid overlay of only one leg at a time:

.. _3now: https://www.rcsb.org/structure/3NOW

.. figure:: images/3now_vs_colabfold.jpg

    Experimental structure in green, `ColabFold`_ model in purple.

.. _ColabFold: https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb

As described above, the sensible approach here is to break the predicted model into two 
domains for molecular replacement. For this tutorial, I've done that using 
`phenix.process_predicted_model`__, which does three things:

* Converts the pLDDT values in the model to estimated B-factors,
* Breaks the model into domains based on analysis of the PAE matrix, and
* Trims off low-confidence loops and tails.

__ https://phenix-online.org/version_docs/dev-4380/reference/process_predicted_model.html

From there, I've run molecular replacement with `Phaser`__, leading to a very high-confidence
solution (Log-Likelihood Gain 2331, translation function Z-score 47). That's where we're going
to start. 

__ https://scripts.iucr.org/cgi-bin/paper?he5368

First, use the ChimeraX **File/Open** dialog to open the file 3now_phaser.1.pdb. You should see 
something like this:

.. figure:: images/mr_solution.jpg

Now, open the predicted model (3now_662d5_unrelaxed_rank_1_model_2.pdb), then align it to the 
C-terminal region of the MR solution with:

`match #2 to #1:400-900`__

__ cxcmd:match\ #2\ to\ #1:400-900

.. figure:: images/mr_solution_af_overlay.jpg

This second model is going to become our working model. If you haven't already started ISOLDE,
launch it now:

`isolde start`__

__ cxcmd:isolde\ start

\... and tell it to focus on this model either by using the "Working on:" drop-down menu, or via
the command:

`isolde select #2`__

__ cxcmd:isolde\ select\ #2

Before we actually load the MTZ file, there's one very important thing we need to do. Since this 
is an AlphaFold model the B-factor column is currently filled with pLDDT values, which are no good 
for crystallographic calculations. Luckily, as described in the supplementary information of 
`the original RoseTTAFold manuscript`__, there's a fairly strong relationship between pLDDT and 
experimental B-factors, allowing us to estimate the latter. On the top ChimeraX menu, go to 
**ISOLDE/Prep/Convert pLDDTs to B-factors** to make this conversion.

__ https://doi.org/10.1126/science.abj8754

Once that's done, we can go ahead and open the MTZ file provided by Phaser. The easiest way to 
do that is via the "Add map(s) to working model" widget on the ISOLDE GUI:

.. figure:: images/add_maps_widget.png

Click the "From crystallographic dataset" button and choose 3now_phaser.1.mtz.

*(NOTE: Because the Phaser MTZ file does not contain free flags, ISOLDE will choose a new set at this 
point.)*

Your display should now look something like this:

.. figure:: images/mtz_loaded.jpg

The R-factors calculated from the experimental structure factors should appear briefly in the bottom
right corner (don't worry if you don't see them yet). They're currently hovering around 0.6 - telling 
us that (due to the n-terminal portion of the model being way out of position) the maps being generated
by ISOLDE itself are currently far from useful. That's OK, since the MTZ file also contains a "static" 
map pre-calculated by Phaser from the MR solution (well, actually it contains three, but we're only 
interested in one). Open the "Precalculated Crystallographic Map Settings" widget on ISOLDE's General tab
to take a look:

.. figure:: images/precalc_maps_widget.png

The three maps here are, in order:

------------------ ---------------------------------------------------------
Name               Meaning
------------------ ---------------------------------------------------------
FC, PHIC           Structure factors calculated from the MR model. This 
                   map will *always* look like the starting model, and is
                   not useful for our purposes.

FWT, PHWT          The 2mFo-DFc structure factors (Two times the observed 
                   amplitudes minus the calculated amplitudes combined with 
                   phases from the calculated map; the "m" and "D" are 
                   weighting factors). This is the map we want to work 
                   with.

DELFWT, PHDELWT    The mFo-DFc map. Not enormously useful to us in this 
                   case.
------------------ ---------------------------------------------------------

Let's delete the maps we don't need:

`close #2.1.1.5|#2.1.1.7`__

__ cxcmd:close\ #2.1.1.5|#2.1.1.7

\... and enable the remaining map for fitting by checking the MDFF checkbox. 
This will pop up a warning message:

.. figure:: images/precalc_map_mdff_warning.png

In general you should pay close attention to this message: if for any reason 
you choose to do the majority of your fitting with a precalculated map rather 
than ISOLDE's own live maps, it is up to you to make sure that map was 
calculated without free reflections included. In this case, though, we'll 
only be using this map fairly briefly for initial fitting, so there's very
little risk here. Go ahead and click OK.

Personally, I also prefer to switch from mesh to transparent surface visualisation
here, but that's up to you.

Now, collapse the precalculated maps widget and open the "Dynamic Crystallographic
Maps" widget above it. It should look like this:

.. figure:: images/dynamic_maps_widget.png

Since these aren't much use to us just yet, we should disable MDFF here and 
hide all the maps to reduce visual confusion:

.. figure:: images/dynamic_maps_hidden.png

Now you can collapse that widget if you like. We're almost ready to start now - 
we just need to add hydrogens:

`addh #2`__

__ cxcmd:addh\ #2

\... and add some restraints to reinforce the model during this initial fitting
(without these, models this far out of density tend to get pretty mangled). For
the restraint reference, let's open a copy of the reference model
(**File/Open**, and choose 3now_662d5_unrelaxed_rank_1_model_2.pdb).

Switch to ISOLDE's Restraints tab, and expand the "Reference Models" widget.
Choose the model you just opened in the "Reference model:" drop-down menu. The
widget should now look like this:

.. figure:: images/reference_model_no_PAE.png

Click the "Load PAE matrix" button, choose the file
3now_662d5_unrelaxed_rank_1_model_2_scores.json, and click OK. In the "Assign
restraints" window, check the checkboxes under Distances and Torsions. The
widget should now look like this:

.. figure:: images/reference_model_with_PAE.png

Feel free to expand the Options and play around with the settings, but when
you're done, make sure to leave it looking like this before continuing:

.. figure:: images/reference_model_options_widget.png

Click the Apply button to add distance and torsion restraints to your working
model. The main window should now look like this:

.. figure:: images/restraints_applied.jpg

One last thing: with the default mask radius of 4 Angstroms (applied when starting
a simulation) the N-terminal portion of the map is too far from our working model 
to be shown. Increase it to 9 Angstroms (On ISOLDE's General tab).

**(NOTE: very large mask radius values can significantly impact on simulation 
performance, particularly when working with live maps. Don't forget to set it back
to a smaller value once this bulk fitting task is done.)**

Now, go ahead and start a simulation. Select the model and press the play button 
on the left of ISOLDE's ribbon menu, or:

`sel #2; isolde sim start sel`__

__ cxcmd:sel\ #2;isolde\ sim\ start\ sel

Click the pause button once the simulation starts. Now, for big bulk-movement tasks 
like this the default all-atom view can make it a bit hard to see the forest through 
the trees:

.. figure:: images/sim_started.jpg

\... so I'd recommend reducing to a C-alpha trace by hiding all other atoms:

`hide ~@CA`__

__ cxcmd:hide\ ~@CA

.. figure:: images/ca_trace.jpg

Now, select the first 300 residues:

`sel #2:1-300`__

__ cxcmd:sel\ #2:1-300

\... and choose the "Tug selection" mouse mode on ISOLDE's ribbon menu:

.. figure:: images/tug_selection_mode.png

Your *right-click-and-drag* tugging actions will now be distributed across all selected 
atoms. Resume the simulation, and use this to help the N-terminal domain into place.
It's generally best to use a series of short tugs with breaks in between rather than 
one aggressive action - your goal is to "help" the model into place, not force-fit it.
It should only take 3-4 such short tugs for it to fall into place.

Before:

.. figure:: images/tug_selection_start.jpg

After:

.. figure:: images/tug_selection_end.jpg

Once your simulation looks like the latter image, hit the green stop button. Take a look 
at the bottom right of the ChimeraX window - the R-factors, previously hovering around 
0.6, should now be down to the vicinity of 0.4. Don't worry if the R-work is above R-free
(hinting at overfitting) at this stage - we're about to drop the Phaser map and switch 
to ISOLDE's live potential - which strictly excludes the free set - for all further rebuilding,
which will rapidly wipe out any model bias.

.. figure:: images/r_factors_after_bulk_fit.png

So let's go ahead and make that switch-over. Disable MDFF on Phaser's FWT,PHWT map, and hide 
it:

.. figure:: images/disable_precalc_mdff.png

Actually, if you prefer you can just close it entirely - it's not needed from here on:

`close #2.1.1.6`__

__ cxcmd:close\ #2.1.1.6

\... display all dynamic maps except the MDFF potential, and enable MDFF on the latter.

.. figure:: images/enable_dynamic_mdff.png

If you wanted to, you *could* display the MDFF potential as well, but in general it's not 
that useful. It differs only from the 2mFo-DFc map you see in that it excludes the free
reflections - the displayed maps should be slightly more informative since they include all 
reflections.

At this stage we can also close the original MR model:

`close #1`__

__ cxcmd:close\ #1

Before we start our next simulation, if you're like me you'll probably be finding that cobweb
of distance restraints a bit hard to see through. To cut down on the clutter, we can tell 
ISOLDE to only show those that are stretched away from their target values. You can do this 
on ISOLDE's "Manage/Release Adaptive Restraints" widget on the Restraints tab - just drag 
the "Display threshold" slider to the right:

.. figure:: images/adjust_restraint_display_threshold.png

This widget also contains some useful tools for selectively releasing distance and torsion 
restraints - we'll be needing these in a bit, because there are various sites where the 
AlphaFold model disagrees with the experimental density. Some of these disagreements are 
probably simply due to the slightly different overall conformation and the lack of crystal 
contacts - but some do look like true errors on AlphaFold's part. Before we get into looking
through all that in detail, let's first run a quick simulation in the new map:

`sel #2; isolde sim start sel`__

__ cxcmd:sel\ #2;isolde\ sim\ start\ sel






