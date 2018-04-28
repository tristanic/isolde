Adjusting map weightings
========================

.. contents::
    :local:

Overview
--------

Weighting of maps (that is, deciding how hard they "pull" on atoms compared to
the strength of geometric restraints) can be a challenging problem. This is
particularly true of cryo-EM maps: whereas crystallographic maps are calculated
in a consistent unit system (electrons per cubic Angstrom), the units of
cryo-EM maps are somewhat arbitrary and can vary over orders of magnitude from
map to map. Furthermore, there is not yet an agreed standard method for their
normalisation. At present, ISOLDE makes no attempt to optimise map weighting -
this choice is left up to the user.

That being said, the default weighting of
:math:`1000\ kJ\ mol^{-1}\ (\text{map density unit})^{-1}\ \AA^3` for standard
maps and :math:`500\ kJ\ mol^{-1}\ (\text{map density unit})^{-1}\ \AA^3`
for difference maps tends to work well for most crystallographic maps. In
general, you should expect to reduce this weighting as the resolution decreases.

Choosing your map weighting
---------------------------

To choose a suitable weighting for your particular map, I suggest the following
protocol:

    1. Navigate to a well-resolved region, and adjust your map contours to give
       give good visual discrimination. Note the density value at this contour
       level (the first number that appears in the bottom left corner of the
       ChimeraX main window as you scroll).
    2. Set your initial guess for the map weighting (click the
       *Show map settings dialogue* button on the *Sim settings* tab) such that
       (weighting) * (optimal contour level) ≈ 1.0 (halve this if you have two
       overlapping maps - e.g. a standard and a sharpened map - both pulling on
       your model).
    3. Make a modestly-sized selection in the vicinity (e.g. one well-resolved
       helix/loop), and start an interactive simulation.
    4. Using the right mouse button, drag an atom ~2-3 \AA out of the map, and
       observe what happens. If it doesn't readily fall back into place, your
       weighting is too low. If the mouse tugging is too weak to pull the atom
       out of the map, your weighting is too high.
    5. In addition to the above, pay attention to the general behaviour of the
       mobile atoms. The appearance of Ramachandran outliers in the middle of a
       well-resolved helix, for example, or flexible sidechains "folding back"
       on themselves to fall into the density is a very good sign that your
       weighting is too high.
