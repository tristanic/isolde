"" 
.. _alphafold_reference_tutorial:

Improving Legacy Models Using AlphaFold
=======================================

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

Tutorial: Quickly improving 2rd0 with AlphaFold reference restraints
--------------------------------------------------------------------

.. toctree::
    :maxdepth: 2

**NOTE: This tutorial assumes that you are already familiar with the basic
operation of ISOLDE. If this is your first time using it, it is highly
recommended that you first work through at least the** :ref:`isolde_intro_tutorial`
**tutorial before starting this one.**

It is no big secret that models built into lower-resolution datasets have 
historically been fraught with errors, but it can often be surprising to the 
uninitiated just how big these errors can get. In the introductory tutorial we 
looked at rebuilding a small 3.0 angstrom model, 3io0, in ISOLDE using only the 
information available to us from the experimental data. While that case was fairly
straightforward, for larger structures, lower-resolution and/or poorer-quality 
starting models this approach still becomes quite time-consuming and occasionally 
frustrating. Luckily, with the advent of AlphaFold and its high-quality models of 
most well-structured regions of proteins, we now have a much better way to quickly
and automatically fix most of the problems in legacy models, reducing the interactive
component to the (typically few) sites where the AlphaFold model disagrees with the 
data, and/or the reference restraints need a little help steering the experimental 
model to the new conformation.

As our example, we're going to take a look at `2rd0`_, a 3.05 angstrom structure of the 
human p110alpha/p85alpha complex. Published back in 2007 as the first 

.. _2rd0: https://www.rcsb.org/structure/2RD0

