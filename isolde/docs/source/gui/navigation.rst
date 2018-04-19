Visualisation and Navigation
============================

.. contents::
    :local:

*(General note: if you open a large model in ChimeraX it will often
automatically switch to its high-quality rendering mode with full shadows and
ambient occlusion. While beautiful, this is sadly too slow to work with ISOLDE.
You can switch back to simple lighting by clicking the yellow spotlight button
in ChimeraX's top panel, or entering "lighting simple" in the ChimeraX command
line)*

General layout
--------------

After loading a model and map (or clicking "Load demo"), your ChimeraX display
should look something like the screenshot below (you may need to use the
ChimeraX atom display controls to set the atom display styles to your liking
first).

.. figure:: images/loaded_model.jpg
    :alt: Model/Map display in ChimeraX window

    A typical ISOLDE scene

First, let's talk about what's changed in the model itself. Perhaps most
immediately obvious is the change in the appearance of the cartoon: it's now
much thinner than you're probably used to seeing. This is to ensure that it
doesn't get in the way of seeing the atoms themselves, while still providing
valuable information about overall topology and secondary structure.

Next (if you're working with a crystal structure) you might notice that your
model has been joined by one or more darker copies of itself. These are of
course the symmetry contacts in the crystal lattice. The symmetry atoms are
non-interactive "ghosts" - while they will update instantly when the "real"
atoms change, you cannot select or move them. Hovering your mouse over one,
however, will give you a popup telling you its name and symmetry operator:

.. figure:: images/symmetry_tooltip.png
    :alt: Symmetry tooltip

    Symmetry atoms know who they are

You'll also note that you can no longer see all of the atoms (nor all of the
map). By default atom display is restricted to residues approaching within 15Å
of the central pivot point, while the map is restricted to a 12Å sphere. Other
display options suited to isolating issues in low-resolution maps will be
discussed  below.

Zooming and Panning
-------------------

*(NOTE: Some ChimeraX functions may change the behaviour of the centre of
rotation to a mode incompatible with the behaviour described below. If you
find things not behaving as they should (e.g. the pivot indicator no longer
remains in the centre of the screen), type "cofr center showpivot true" in the
ChimeraX command line or just click the ISOLDE Spotlight Mode described below.)*

If you've spent some time using ChimeraX before, you've probably already tried
to zoom in using the scroll wheel. That won't work in ISOLDE: since model
building requires regular adjustment of map contours, the scroll wheel is
co-opted to perform that all-important function. A special zoom mode (designed
to bring you inside the model while fading out background details) has instead
been mapped to **shift-right-click-and-drag**. Panning (that is, translating
the display up-down and left-right) is the usual **middle-click-and-drag**.

Adjusting the map contours
--------------------------
