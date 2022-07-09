.. _alphafold_multimer_cryoem_tutorial:

Building Into CryoEM Maps with AlphaFold Models
===============================================

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


Tutorial: Fitting a new cryoEM map starting from an AlphaFold multimer prediction
---------------------------------------------------------------------------------

*Level: Moderate*

.. toctree::
    :maxdepth: 2

**This tutorial will write a number of files to your current working directory. To 
avoid clutter and confusion it is best if you change to a fresh directory first, using
the command:** 

cd browse 

**Due to the way ChimeraX tutorials are implemented you will 
need to enter this command in the ChimeraX command line yourself.**


Historically, the typical approach to building an atomic model into a new crystallographic
or cryoEM map was "bottom-up": starting either from nothing or some partial model representing
the rigid, conserved core, the chain(s) would be traced through the density by some combination
of automated algorigthms and manual building. With the advent of largely-reliable machine
learning-based structure predictions heralded by the release of AlphaFold, it is now usually
both faster and more reliable to instead take a "top-down" approach, starting from (the) 
AlphaFold prediction(s) and refitting them into the new experimental map.

For this tutorial, we are going to largely recapitulate 7drt_, the 2.2 Angstrom structure of 
the human Wntless:Wnt3a complex. `Published in 2020`_ and released in 2021, this is the first 
representative of this complex in the wwPDB and was not included in the training set for AlphaFold2
or the later multimer implementation, making it an excellent test case to explore the strengths 
and weaknesses of this approach. 

.. _7drt: https://www.rcsb.org/structure/7DRT
.. _Published in 2020: https://dx.doi.org/10.1038/s41467-021-24731-3 

Let's start by taking a look at the original model and map:

`open 7drt; open 30827 from emdb; color bychain`__

__ cxcmd:open\ 7drt;open\ 30827\ from\ emdb;color\ bychain

\... and make it a little more pretty:

`volume #2 rmsLevel 8 step 1;lighting soft`__

__ cxcmd:volume\ #2\ rmsLevel\ 8\ step\ 1;lighting\ soft

.. figure:: images/7drt_original_overview.jpg

The first thing you might note when looking at this is that, as in many cryoEM reconstructions, that
"2.2 Angstrom" resolution is a nominal value relating to the best-resolved parts of the map (in this
case the transmembrane domain), and that other regions have significantly lower resolution - revealed
as noisy and/or absent density at this contour level. The resolution around the top of the Wnt3a chain 
is particularly poor - in sites like this we'll have to rely heavily on prior knowledge.

Preparing the predicted model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Our first step in building into this map based on a predicted model is, of course, to actually generate 
the model. There are now many tools you can use for this, most of which make use of Google Colab so you
don't need to have an AlphaFold installation on your own computer. Examples include 
`the official DeepMind version`_, the `ColabFold version`_ using a faster multiple sequence alignment
implementation, and the `ChimeraX version`_ directly launchable from within ChimeraX.

.. _the official deepMind version: https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb

.. _ColabFold version: https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb

.. _ChimeraX version: cxcmd:ui\ tool\ show\ AlphaFold

In my experience each of these yields quite similar results, so your choice of which to use is mostly a 
matter of preference. Except for very small "toy" examples they're still too slow to run in the context 
of an interactive tutorial, so here we'll be working with the stored results from a ColabFold run (which
took about 1.5 hours to generate 5 independent predictions of this complex).

Before we get started properly on this, we need to consider a few caveats. While these predictions are often 
startlingly good, for various reasons they are rarely perfect. The first and most obvious of these reasons is 
that the prediction only gives one conformation - even if this correctly represents a natural state of the protein,
that is not necessarily the same as *your* state. Other potential issues include:

* AlphaFold2 was trained on essentially the entire protein cohort in the wwPDB, with very little filtering 
  for resolution or model quality. As a result, it sometimes makes similar errors in sidechain geometry to those 
  commonly seen in human-generated experimental models.
* In places where there is too little information to control things (generally because a site is intrinsically 
  disordered) the local geometry can rapidly degrade into nonsense. It's important to make sure we discount 
  sites like this.
* For multimers like this where there is no pre-existing structure in the training set, it has to rely entirely on sequence coevolution
  combined with whatever it has learned about what "belongs together" to guide the fold. Importantly, due to the 
  design of the multimer algorithm it **only works for complexes where all chains are from the same organism**,
  and does best on tightly-conserved "exclusive" complexes (those with few or no non-interacting paralogs). Don't 
  expect a sensible result if you try to predict an antibody:antigen complex, for example.
* AlphaFold2 knows nothing of ligands or post-translational modifications (or of physics, for that matter). Nevertheless,
  it will often recapitulate geometry with space for the missing ligands.

Anyway, let's go ahead and get started. First, we need to copy the predicted models to our tutorial directory. If you 
haven't already changed to a clean directory with "cd browse", do that now first. Then, run the following command 
*(Due to the details of ChimeraX help file implementation, you'll have to type or copy-paste this yourself)*:

isolde tut prep alphafold_cryoem_multimer

This will create a new directory '7drt' in your working directory with four files in it (two models and two JSON files 
containing their Predicted Aligned Error (PAE) matrices). A real run of ColabFold typically creates 5 models and returns 
a few other useful files, but these are omitted here to save space. Anyway, go ahead and open the top-ranked model 
("7drt_20b0f_unrelaxed_rank_1_model_4.pdb") via the ChimeraX File/Open dialog. 

Docking into the map
~~~~~~~~~~~~~~~~~~~~

Before doing anything that involves moving your model, it's best to reduce the lighting effects back to the simple mode - 
ambient occlusion is beautiful, but expensive to calculate!

`lighting simple`__

__ cxcmd:lighting\ simple

At present ChimeraX does not provide an algorithm that will dock a model into a map starting from a random position, so you 
will need to do the initial work yourself. Of course in this case we could cheat by aligning to the existing experimental 
structure with

`match #3 to #1`__

__ cxcmd:match\ #3\ to\ #1

\... but the more general solution is to get most of the way with the "Move Model" right mouse mode, then use *FitMap* to 
optimise the fit. Go to the "Right Mouse" tab on the ChimeraX ribbon menu and click the "Move model" button (fifth from the 
right), or use the command:

`ui mousemode right "translate selected models"`__

__ cxcmd:ui\ mousemode\ right\ "translate\ selected\ models"

Then select the AlphaFold model:

`sel #3`__

__ cxcmd:sel\ #3

and to make things fair, hide the experimental model:

`hide #1 model`__

__ cxcmd:hide\ #1\ model

Now, *right-click-and-drag* will move the model across the screen. Pressing shift while dragging (**AFTER** clicking and holding the 
right mouse button) changes this to rotation. Left-click-and-drag rotates the display, and scroll zooms. Use these actions to bring 
your model approximately into line with the density. Once you think you're close enough, do:

`fit #3 in #2`__

__ cxcmd:fit\ #3\ in\ #2

If the result isn't to your liking, shift the model a bit and try again. You'll probably notice quite quickly that it's impossible to 
get a satisfactory fit to the transmembrane domain, because the helices are shifted substantially relative to each other. Instead, focus
on optimising the fit of the big beta-sandwich domain of Wntless:

.. figure:: images/af_model_1_fitted.jpg


You can limit *fitMap* to only consider this domain with:

`fit #3/C:50-230 in #2`__

__ cxcmd:fit\ #3/C:50-230\ in\ #2

Don't worry that the fit isn't perfect - after all, that's what we're here to fix! But credit where it's due - 
considering that when AlphaFold2 was trained there was no experimental structure for this complex (or any homologue)
this model is *amazingly* good. Anyway, if at least part of the model is well-fitted,
that makes everything following that much easier.

At this point, it's finally time to initialise this model for ISOLDE. First, let's close the experimental model:

`close #1`__

__ cxcmd:close\ #1

If you haven't already launched ISOLDE, start it now:

`isolde start`__

__ cxcmd:isolde\ start

Then, click the "Working on:" button at the top left of the ISOLDE window, or alternatively use the command:

`isolde select #3`__

__ cxcmd:isolde\ select\ #3

A dialog box will appear noting that the model looks like it should have disulfide bonds that aren't specified 
in the metadata and asking if you wish to create them. Choose "Yes". 

**NOTE: when first prepping a model, ISOLDE temporarily removes it from the ChimeraX session before returning it
under the control of a higher-level "manager" with the first available model ID. In this case, the working model 
now becomes model #1.**

Before we associate the map with the model and start rebuilding, it would be a good idea to take a quick look 
at its quality. For clarity, let's temporarily hide the map:

`hide #2 model`__

__cxcmd:hide\ #2\ model

Now, go to ISOLDE's "Validate" tab and open the Clashes widget.

.. figure:: images/model_1_clashes.png

Hmmm - that's not a good sign. A carbon-carbon overlap of 2.5 Angstroms means the atoms are occupying almost 
the same position, which could be a real headache to untangle. Click on the top entry in the table to take a 
closer look. 

.. figure:: images/bad_clash.jpg

Ugh. That's not good. This looks like an instance of one of the known deficiencies of AlphaFold: when two stretches
(a) are distant in sequence space and (b) have no mutual information controlling their relative positions, in key
stages of the algorithm they don't "know" about each other and can end up overlapping like this.

While a situation like this *is* recoverable with judicious use of the "isolde ignore" command (excluding one of the 
offending stretches from simulation while clearing the other), let's see if we can avoid the issue by using the second 
model instead. Go ahead and open that ("7drt_20b0f_unrelaxed_rank_2_model_5.pdb") with the File/Open dialog, and align
it to the working model:

`match #3 to #1`__

__ cxcmd:match\ #3\ to\ #1

You might find it easier to compare if you temporarily hide all atoms:

`hide #1`__

__ cxcmd:hide\ #1

.. figure:: images/model_2_no_clash.jpg
  








 