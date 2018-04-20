Basic simulation controls
=========================

.. contents::
    :local:

Starting a simulation
---------------------

Assuming your model meets the reuqirements (that is, all residues are complete
with hydrogens and known to the MD forcefield), then starting a simulation is as
simple as making a selection of atoms in the main ChimeraX window, then clicking
the play button in the bottom left of the ISOLDE panel (see below). Your
selection may be anything from a single atom to all atoms in the model, but for
the sake of performance should be as small as necessary to achieve the desired
goal. Of course, the larger the selection the slower the simulation.

The selection size you can make while still achieving reasonably interactive
performance is primarily dependent on your GPU (which does essentially all of
the hard work of running the simulation). As a general reference, my "big"
laptop (a gaming model with Intel Core i7 CPU and a NVidia GTX1070 GPU)
manages about 15-20 coordinate updates per second when simulating all 229
residues of the demonstration model, and maintains reasonably interactive
performance for up to a few thousand residues. Note that *graphics* performance
(that is, responsiveness of the display when rotating, translating, zooming
etc.) is only weakly affected by simulation size, and remains above 30 fps
under most circumstances. My (much) smaller MacBook Air still manages
surprisingly well. Not that I'd actually suggest seriously using such a small
machine for an entire job, but it gives usably interactive performance for a
few dozen residues at a time (sufficient for tweaking rotamers, loops etc.).

The maximum sized model that will successfully run **non** interactively is a
limit I have yet to find, but will typically be limited by the amount of RAM on
your GPU. I can tell you that simulations of the 3,792-residue 3ja8 with its
3.8Å map run successfully on both of the above machines. While in most cases
simulation of the entire model will be far too slow to use interactively, it is
nevertheless often useful to run it for a few minutes prior to running any
interactive simulations on smaller selections, to ensure any bad clashes and
other very-high-energy states are relaxed out.

Before we go ahead and start a simulation, let's talk briefly about what happens
to your initial selection when you press that play button. That's controlled by
this dialogue on the *Sim settings* tab:

.. figure:: images/sim_selection_settings.png
    :alt: Simulation selection expansion settings

    Control how the selection will be expanded for simulation

The first spinbox controls how the selection will be extended *along* the
chain(s), while the second controls how the resulting selection will be expanded
*outwards*. The steps taken are:

    1. The initial selection is expanded to complete residues (for any residue
       in which *one* atom was selected, *all* are selected);
    2. For all residues in linear chains (i.e. protein or nucleic acid), each
       contiguous selection is extended forward and backward by the number of
       residues given in the first spinbox, stopping at chain breaks.
    3. Any residue for which any atom comes within the distance given in the
       second spinbox from any atom in the selection defined by (2) is added to
       the selection.
    4. The selection defined by (3) contains all the residues that will be
       mobile in the simulation. In addition, a shell of rigidly-fixed atoms
       is added surrounding the mobile selection, to make sure all mobile atoms
       maintain their physical context in the wider model.
    5. If a distance restraint involves a mobile atom and a second atom which
       *isn't* in the mobile or fixed selections, the residue containing the
       second atom is added to the fixed selection.

The default values are fine for most purposes, but feel free to play. Just keep
in mind that *very* small simulations with only a few mobile residues
surrounded by fixed atoms tend to be somewhat unstable, since they often have
no freedom to relax bad (high-energy) interactions.

General simulation controls
---------------------------

Anyway, go ahead and click the play button, and let's go through the rest of the
buttons in this toolbar:

.. figure:: images/simulation_controls.png
    :alt: Simulation control toolbar

    General simulation controls

    +----------------+---------------------------------------------------------+
    | |play| |pause| | If no simulation is running, starts one. If a simulation|
    |                | *is* running, controls pause/resume.                    |
    +----------------+---------------------------------------------------------+
    | Min/Equil      | Switch the simulation between energy minimisation and   |
    |                | equilibration modes. (*NOTE: due to the way minimsation |
    |                | works under the hood, minimisation is not very          |
    |                | interactive. A quite effective "interactive             |
    |                | minimisation" can be achieved by simply setting the     |
    |                | temperature to 0 Kelvin in equilibration mode*)         |
    +----------------+---------------------------------------------------------+
    | |temp|         | Set the simulation temperature in Kelvin. The default   |
    |                | 100K provides enough thermal motion to help "jiggle"    |
    |                | atoms into more favourable conformations, without moving|
    |                | so fast as to be hard to control. The temperature can   |
    |                | be adjusted to any value between 0 and 500K.            |
    +----------------+---------------------------------------------------------+
    | |cpg|          | Set a "checkpoint" for the current simulation.          |
    |                | A checkpoint is simply a snapshot of the current        |
    |                | simulation, including the atomic coordinates and the    |
    |                | states of all interactive restraints. You should        |
    |                | consider saving a checkpoint whenever you're about to   |
    |                | embark on a tricky or tentative manipulation of your    |
    |                | model. Note that the checkpoint becomes invalid once    |
    |                | the simulation is stopped.                              |
    +----------------+---------------------------------------------------------+
    | |cpr|          | Return to the checkpoint you last saved, discarding all |
    |                | subsequent changes. If you have not saved a checkpoint  |
    |                | yet, this will revert you to the start of the           |
    |                | simulation.                                             |
    +----------------+---------------------------------------------------------+
    | |stopg|        | Stop the simulation, keeping the current state.         |
    +----------------+---------------------------------------------------------+
    | |cps|          | Stop the simulation and revert to the last checkpoint.  |
    +----------------+---------------------------------------------------------+
    | |stopr|        | Abort the simulation, discarding all chainges (this will|
    |                | raise an "are you sure?" warning).                      |
    +----------------+---------------------------------------------------------+

.. |play| image:: ../images/play_icon.png
.. |pause| image:: ../images/pause_icon.png
.. |temp| image:: ../images/thermometer.png
.. |cpg| image:: ../images/checkpoint_green.png
.. |cpr| image:: ../images/checkpoint_red.png
.. |stopg| image:: ../images/stop_sign_green.png
.. |cps| image:: ../images/checkpoint_stop_red.png
.. |stopr| image:: ../images/stop_sign_red.png

Basic interaction
-----------------

With your simulation running, **right-click-and-drag** on any heavy atom and
you should see it follow your mouse cursor, with a 3D green arrow showing the
direction of pull:

.. figure:: images/tugging.png
    :alt: Tugging on a single atom

    Interactive tugging forces are shown like this green arrow

This mouse mode is set every time the simulation starts. You can still use the
buttons on the ChimeraX mouse modes toolbar to make the right mouse button do
other things (just avoid anything that actually edits the model while your
simulation is running, otherwise Bad Things\ :sup:`TM` will happen), and come
back to tugging using the ISOLDE tugging mouse modes toolbar (at the bottom left
of the ISOLDE panel):

.. figure:: images/tugging_toolbar.png
    :alt: ISOLDE's tugging modes toolbar

    Available tugging modes (expect this to grow over time)

    +-------------+------------------------------------------------------------+
    | |tuga|      | **right-click-and-drag** tugs on a single non-hydrogen     |
    |             | atom.                                                      |
    +-------------+------------------------------------------------------------+
    | |tugr|      | **right-click-and-drag** tugs on all heavy atoms in a      |
    |             | residue. The total applied force is the same as when       |
    |             | tugging on a single atom, but is distributed across all    |
    |             | atoms in proportion to their mass.                         |
    +-------------+------------------------------------------------------------+

.. |tuga| image:: ../images/tug_atom.png
.. |tugr| image:: ../images/tug_residue.png
