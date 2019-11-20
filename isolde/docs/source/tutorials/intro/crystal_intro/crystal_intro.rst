.. _isolde_intro_tutorial:

Introduction to crystallographic model rebuilding in ISOLDE
===========================================================


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

Tutorial: Diagnosing and rebuilding errors in a 3Å crystal structure
--------------------------------------------------------------------

.. toctree::
    :maxdepth: 2

Building models into low-resolution crystallographic density maps has
traditionally been a difficult and often frustrating task. To illustrate why,
let's look at what happens to the density as resolution decreases. You can also
use this as an opportunity to familiarise yourself with the controls.

Let's start by looking at a small, high-resolution (1.1Å) model, 1a0m. The
Clipper plugin allows you to fetch this directly from the wwPDB along with its
structure factors using the  ChimeraX command line:

`open 1a0m structureFactors true`__

__ cxcmd:open\ 1a0m\ structureFactors\ true

After running this command, check the ChimeraX log. You should see a warning
message similar to that below:

.. figure:: images/1a0m_warning.png

In crystallographic model building it is standard practice to leave a random
selection of measured reflections (typically 5% or 2000 reflections, whichever
is higher) as a guard against overfitting (see `Brunger, 1992`__). While it is
now standard practice for structure factors deposited to the wwPDB to include
the free reflections used by the original modellers, this was not always the
case. When no free flags are found a fresh set will be created for you - in
such cases you should follow the advice in this warning and immediately save
a fresh MTZ file incorporating these:

__ https://www.nature.com/articles/355472a0

`save ~/Desktop/1a0m.mtz #1`__

__ cxcmd:save\ ~/Desktop/1a0m.mtz\ #1

*(OPTIONAL: personally I prefer to work with a white background. If you're
the same, you can change it using ChimeraX's icons in the "Home" tab at top,
or with the following command)*

`set bgColor white`__

__ cxcmd:set\ bgColor\ white

Now, zoom in by **scrolling** the mouse, and let's go ahead and look at the
model you just opened.

.. figure:: images/1a0m_sample.jpg

As you can see, at this resolution (almost) every non-hydrogen atom has its own
distinct, clearly-differentiated "blob" of density - that is, the map shows you
precisely where each individual atom belongs, and any errors in the model tend
to be quite clear (and have equally clear solutions).

In this particular case, there **are** no real errors in the model - as one
would hope given the resolution! While ultra-high-resolution models like this
aren't really what ISOLDE was designed for, this one makes for a nice "toy" to
play with while learning the controls.

Getting around
~~~~~~~~~~~~~~

While zoomed in, you'll notice a small set of axes displayed in the middle of
the screen:

.. figure:: images/pivot_indicator.jpg

This is the pivot indicator marking the centre of rotation, and is useful as a
reference point for certain tasks. The colours represent the three primary axes:
red=x, green=y, blue=z.

When you zoom, the distance between the near and far clipping planes
automatically grows and shrinks to provide depth cueing and reduce distracting
clutter. If you need to adjust this, you can do so with **shift-scroll**.

Now, scroll out to a comfortable distance, and try panning. You have two options
here: if you are using a wired mouse, just **click and drag with the middle
mouse button**. For touchpad users without a middle-button equivalent, use
**shift-left-click-and-drag**. To move the centre of rotation in the Z direction
(that is, towards or away from you) use **ctrl-middle-click-and-drag**.

You'll hopefully notice a few things:

1. Wherever you go, you're surrounded by a sphere of density and displayed
   atoms. This is the default navigation mode in ISOLDE. The radius of the
   sphere can be controlled using the ISOLDE GUI (more on that in a bit) or
   using a command like `clipper spotlight radius 15`__

__ cxcmd:clipper\ spotlight\ radius\ 15

2. Zoom out until you have a clear border all around the displayed atoms, then
   pan to the right. As you go, you'll see a single cartoon start to stand out
   from the rest:

   .. figure:: images/separate_asu.jpg

   This is the **actual** model, which is always displayed in this mode -
   everything else is "ghost" drawings of the symmetry atoms (identifiable by
   their darker colour relative to the primary model). While you can't interact
   directly with these, hovering over any symmetry atom will bring up a tool-tip
   giving its identity and symmetry operator. Note: if you ever get so lost that
   you can't find your true model any more, just type `cview #1.3`__
   in the command line, where #1.3 is the model number of the atomic model in
   the model panel:

   .. figure:: images/model_panel.png

__ cxcmd:cview\ #1.3;

The "cview" command is very similar to ChimeraX's "view" command - the only
difference is that "view" fixes the centre of rotation on the centroid of
the specified selection whereas "cview" maintains Clipper's standard mobile
centre of rotation.


Adjusting your maps
~~~~~~~~~~~~~~~~~~~

When you load a model with structure factors, the Clipper plugin generates three
maps: a standard 2mFo-DFc map; a second 2mFo-DFc map with either sharpening or
smoothing applied depending on resolution (maps better than 2.5Å are smoothed,
lower resolutions are sharpened); and a mFo-DFc map. The sharpest of the
2mFo-DFc maps is displayed as a transparent cyan surface, and the other as cyan
wireframe. The sharper map gives the best detail in the well-resolved portions of
the model, while the smoother one is better at revealing connectivity in poorly-
resolved regions where sharp maps are often dominated by noise and difficult to
interpret.

The difference map is shown as green and red wireframes for the
positive and negative contours respectively.

On loading, contour levels are automatically chosen that are generally good for
visualising the majority of the model - but you will of course need to adjust
them from time to time. You can contour a given map using **alt-scroll**, and
choose the map you wish to contour using **ctrl-scroll**.

In many cases (particularly when working in low-resolution maps) you may wish
to focus on some extended selection, where the default "sphere of density" view
becomes impractical. You can expand and mask the map to cover any arbitrary
selection using the **clipper isolate** command. For example:

`clipper isolate #1/A&protein`__

__ cxcmd:clipper\ isolate\ #1/A&protein

... should give you something like this:

.. figure:: images/isolate_A.jpg

More specific selections can be specified using ChimeraX's very powerful
`atomspec`__ syntax.

__ help:user/commands/atomspec.html

You can further tweak the results using various optional arguments. To see the
available options, use `usage clipper isolate`__ and check the Log window.

__ cxcmd:usage\ clipper\ isolate

Once you're done, go back to the scrolling-sphere mode using `clipper spotlight`__.

__ cxcmd:clipper\ spotlight

Running your first simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let's go ahead and start ISOLDE. You can do this from the menu
(**Tools/General/ISOLDE**) or by simply running the command `isolde start`__.

__ cxcmd:isolde\ start

You will notice is that all the C-alpha atoms in the model have suddenly turned
into green spheres:

.. figure:: images/rama_markup.jpg

This is one example of ISOLDE's real-time model validation, indicating each
residue's Ramachandran score (that is, the prior probability of the backbone
"twist"  around that residue). You can read more about this `here`__.

__ ../../gui/getting_started.html

If for any reason you wish to turn off display of any of ISOLDE's validation
or restraint markup, you can do so via the Models panel:

.. figure:: images/model_panel_with_validation.png

Just expand the drop-down list under the atomic model, and click the relevant
checkbox under the eye icon. I would recommend against this in most cases,
though: the information provided by the validation markup becomes very useful in
guiding rebuilding, and there is nothing quite as frustrating as fighting
against what turns out to be a restraint that you've hidden! If you wish, you
can tell ISOLDE to show Ramachandran C-alpha markup for only non-favoured
residues using the command:

`rama showFavored false`__

__ cxcmd:rama\ showFavored\ false

... but personally, I prefer to keep them present. If you're the same, put them
back with:

`rama showFavored true`__

__ cxcmd:rama\ showFavored\ true

Before we actually get things moving, it's worth considering some of the
strengths and pitfalls of molecular dynamics methods. Unlike minimisation using
`traditional Engh and Huber restraints`__, molecular dynamics attempts to
realistically model **all** atomic interactions (including the van der Waals and
electrostatic forces between non-bonded atoms) according to Newton's laws of
motion. In most ways this is extremely useful (as we will see), but the results
can be problematic when it comes to very small ligands such as water molecules
and monoatomic ions. In the absence of any explicit bonds to their neighbours,
they can be prone to "wandering". For this reason, ISOLDE provides a command
allowing you to add explicit restraints to all small ligands in the model:

__ http://journals.iucr.org/a/issues/1991/04/00/li0061/li0061.pdf

`isolde restrain ligands #1`__

__ cxcmd:isolde\ restrain\ ligands\ #1

.. figure::images/restrain_ligands.jpg

**(NOTE: Validation and restraint markups are only drawn for real atoms, not
their symmetry ghosts.)**

By default, polar ligand atoms will receive distance restraints (represented by
a cyan bar) to all other polar atoms within 4Å, and ligand residues with fewer
than four distance restraints will be supplemented with position restraints
(represented by a yellow pin). Note that these restraints are defined by the
starting geometry (that is, they make no attempt to impose ideal geometry on the
model), but are soft enough to allow substantial relaxation. Individual
restraints can also be added, adjusted or removed using tools on ISOLDE's
"Rebuild" tab.

There is one more task we need to do before starting our first simulation:
adding hydrogens. Like most molecular dynamics methods, ISOLDE requires residues
to be chemically complete (with a few exceptions: "bare" N- and C-termini and
sidechain truncations are permitted). ChimeraX's *AddH* command typically does a
great job of adding hydrogens. Just run:

`addh`__

__ cxcmd:addh

To reduce clutter in the display, I prefer to hide the non-polar hydrogens:

`hide HC`__

__ cxcmd:hide\ HC

If you forget to add hydrogens before attempting to run a simulation, it's not the
end of the world. ISOLDE will bring up an error dialog like this:

.. figure:: images/sim_error.png

If you receive this error **after** adding hydrogens, stop and carefully inspect
the offending residue. Possible reasons and their suggested remedies include:

* The residue in question is an unusual one for which ISOLDE lacks force field
  parameters. Parameterising new residues is beyond the scope of this tutorial.
  For now, you may choose to ignore the residue (take care: it will remain
  fixed in space, and surrounding atoms will pass through it as if it isn't
  there), or delete it using ChimeraX's **delete** command.

* The residue is missing some of its heavy atoms (or, for an amino acid
  residue, truncated in a way that ISOLDE doesn't support). Truncated
  non-amino-acid residues are currently not supported. Amino acids other than
  proline, threonine, valine or isoleucine may be truncated back to CA, CB, CG
  or CD. If you run into trouble with an unsupported amino acid truncation,
  ChimeraX's top "Right Mouse" tab includes a handy "Swapaa" right mouse mode
  allowing amino acid mutation by right-clicking on an atom and dragging up and
  down. Simply "mutate" the residue to its original identity to fill in all
  missing atoms.

* AddH has added too many or too few hydrogens. Common scenarios include:

  - a peptide bond nitrogen too close to a metal may not be protonated. In
    that case, try repeating the **addh** command with the added argument
    **metalDist 1.0**.

  - **addh** will often add a proton to a terminal phosphate (e.g. on ADP
    and ATP residues). Just select the offending hydrogen using **ctrl-click**
    and type **delete selAtoms** in the command line.

  - for non-protein, non-nucleic acids **addh** does its best to infer
    chemistry from the ligand geometry, which can of course fail if the
    ligand is poorly modelled. ISOLDE does not yet provide a remedy for this
    situation. You may choose to ignore or delete the residue as above, or
    correct the geometry externally before trying again.

In this case, though, the simple **addh** command "just works" - so let's go
ahead and get the simulation moving. First, select the whole model:

`select #1`__

__ cxcmd:select\ #1

... and either click ISOLDE's blue "play" button (at the bottom left of the
ISOLDE panel) or use the command:

`isolde sim start sel`__

__ cxcmd:isolde\ sim\ start\ sel

**(NOTE: you do not HAVE to select the entire model to start a simulation in
ISOLDE. As long as you have at least one atom selected, it will create a
simulation with sufficient scope to encompass your region of interest surrounded
by a soft "buffer zone" giving it freedom to work. In large models you should
expect to spend most of your time in such "mini-simulations", but - particularly
at low resolutions - you should always run at least one "whole model" simulation
to give any clashes and other very-high-energy states a chance to relax.)**

After a few seconds, you should see your atoms start jiggling and the maps begin
to regularly update. When working with crystallographic structure factors in
ISOLDE, every change to the atomic coordinates, B-factors or occupancies
triggers an automatic recalculation of structure factors, maps and R-factors.
A running summary of the current R-factors is provide in the status bar at
bottom right.

Interacting with the model
~~~~~~~~~~~~~~~~~~~~~~~~~~

Now, try tugging on an atom using the right mouse button.

.. figure:: images/tugging.jpg

*(NOTE: Non-polar hydrogens are not tuggable in ISOLDE. Polar hydrogens may be
tugged, but in order to avoid instability the applied force will be much weaker
than that applied to heavy atoms. While tugging on a hydrogen can be useful for
adjusting the H-bonding of a hydroxyl group or water (for example), to make any
larger changes you should choose a heavier atom to pull.)*

This should result in something like the above. There are a few things to note
here:

* The change in the maps (the large red blob around the sulphur atom, and the
  green blob where it **used** to be) is clearly telling you this is wrong
  (quite apart from the R-factors, which in this case have spiked from about
  0.21 to 0.29);
* Tugging on the sulphur has had a knock-on effect on the adjacent Asn sidechain
  (top) - to re-iterate, in ISOLDE, all atoms feel each other.

  **(IMPORTANT NOTE: during simulations, the model will only feel the map as it
  was at the moment the simulation started. Subsequent changes will not take
  effect until you stop the simulation and start a new one.)**

Unless you've done something particularly drastic, upon releasing the tugging
force the model should snap back into place. In my case, here's the model ~5
seconds later, with Rwork/Rfree back at 0.209/0.220:

.. figure:: images/relaxed.jpg

The right mouse button is automatically mapped to this "tug single atom" mode
every time a simulation is started. At the bottom left of the ISOLDE panel you
will find three buttons allowing you to switch between this and "tug residue" or
"tug selection". These modes are primarily useful at much lower resolutions than
this model.

General model controls
~~~~~~~~~~~~~~~~~~~~~~

You may have noticed that ISOLDE's blue play button has changed into a red pause
symbol. Go ahead and click this to pause the simulation while we explore some of
the other available controls.

First, click on the "Show map settings dialogue" button under ISOLDE's "Sim
settings" tab. Here you will find a drop-down menu listing all maps associated
with the model, and tools to adjust both their visualisation and their weight in
the simulation. You will find that for all but one of the maps the weighting
options are disabled. This is because the maps you **see** are generated using
the complete reflection dataset (including the free reflections) and hence
should not be used to guide the model. When working with structure factors,
ISOLDE creates a special map named "(LIVE) MDFF potential" with the free
reflections excluded - choose this from the drop-down menu.

On loading a new map, ISOLDE automatically chooses a suitable weighting based on
the steepness of the gradients in the voxels immediately around the model atoms.
The chosen weighting is usually quite sensible, but you may find you want to
adjust it at times. This is where you can do that. Click the play button to get
the model moving again, then change the weight to 0.3 and click "Set". You
should see an immediate change in the model's behaviour: the atoms will move
much more "loosely", more red and green will appear in the difference map, and
the R-factors should increase to 0.24-0.25.

Injection of random positional error like this is one valid way to reduce the
incidence of model bias that may arise due to the fact that we are not using the
original free set. Clicking the green "STOP" button at bottom right will stop
the simulation and save the instantaneous coordinates. Now, set the weight back
to 3.0 (don't forget to click the Set button!), `select the model`__ and start
a fresh simulation. You should see the R-factors quickly drop back to the
vicinity of 0.21-0.23. Now, drop the temperature to zero using the spinbox
marked by a thermometer icon (to the right of ISOLDE's play button). Wait for
the R-factors to settle to constant values, and hit the green STOP button. Now
start a new simulation - you should see the R-factors drop a little as the
updated potential takes effect.

__ cxcmd:select\ #1

You may choose to continue with the temperature set to zero if you like, but
personally I always prefer to have at least a little thermal motion happening -
if for no other reason that it makes it easy to tell at a glance when a
simulation is actually running! I'd suggest setting the temperature to 20-30K
before continuing.

Now, let's explore the other simulation control buttons in this row. To the
right of the temperature control you'll see two buttons with a green and red
chequered flag respectively. These are the "checkpoint" controls. Clicking the
green flag will save an instantaneous snapshot of the model as it is right now:
not just the atomic positions, but all currently-active restraints. Go ahead and
click it. Now, do something horrible to your model - just grab an atom and pull
as hard as you can. You should see the map explode into a kaleidoscope of red
and green, and the R-factors shoot up into the 0.4-0.6 range:

.. figure:: images/big_mess.jpg

  Whoops.

Click the red chequered flag to put everything back as it was. This is designed
for the not-uncommon situation where the solution to a particular problem isn't
clear, allowing you to try out various hypotheses without the risk of permanent
damage to the model.

To the right of the chequered flag you'll find three different "STOP" buttons.
The first of these (the green one) you've already used - it stops the simulation
and keeps the current state. The next (a red stop button under a red chequered
flag) discards the current state and instead commits the last saved checkpoint
(reverting to the state before the simulation started if no checkpoint was
saved). The big red stop button discards the simulation entirely (after bringing
up a warning dialogue asking if you're sure), returning the model to the state
it was in at the moment the play button was pressed.

Rebuilding
~~~~~~~~~~

While 1a0m is a useful model for learning to find your way around in ISOLDE,
being very high-resolution and close to error-free means there is not really
anything much to rebuild. So, let's close it and move on to a more real-world
task. Make sure you stop you simulation if you haven't already:

`isolde sim stop`__

__ cxcmd:isolde\ sim\ stop

... then close the model:

`close #1`__

__ cxcmd:close\ #1


The model we'll be looking at now is 3io0_, a 3Å model of a small bacterial
shell protein. This is discussed in more detail (including a video covering the
whole rebuilding process) as a case study `on the ISOLDE website`__.

__ https://isolde.cimr.cam.ac.uk/case-studies/case-study-3io0/

.. _3io0: https://www.rcsb.org/structure/3io0

To load these coordinates you could use "open 3io0 structureFactors true" as we
did for 1a0m, but ISOLDE also keeps a locally cached version with hydrogens
already added and an extraneous water molecule removed. To load it, click `here`__.

__ cxcmd:isolde\ demo\ crystal_intro

First, just have a browse around and compare the general appearance of the
density to what you saw in 1a0m. For example, look at residue ASN187:

`cview /A:187`__

__ cxcmd:cview\ /A:187

.. figure:: images/3io0_asn.jpg

For the sake of comparison, here's an asparagine in a similar environment in
1a0m (residue A9):

.. figure:: images/1a0m_asn.jpg

In the high-resolution model, not only is every non-hydrogen atom clearly
resolved but the density even clearly shows the larger atomic radius of the
oxygen compared to the nitrogen. In 3Å density, on the other hand, this detail
is lost: the whole sidechain is reduced to an essentially featureless blob. The
challenge such density poses given traditional restraint schemes is immense:
with only weak guidance from the density itself and none from the non-bonded
interactions with surroundings, finding the *correct* conformation out of all
the possible *incorrect* conformations is, quite simply, hard. It is not
surprising then that models of this resolution have historically been quite
error-prone, and this one is no exception. In fact, a close look at Asn187 in
context with its surroundings suggests it has been built backwards (with the
oxygen occupying the space of the nitrogen and vice versa - look at the
interactions with the peptide bond of Thr303 and the sidechain of Gln199 to see
why).

Anyway, now is a good time to have a look at ISOLDE's "Validate" tab. First,
open up the Ramachandran plot and take a look.

.. figure:: images/rama_plot.png

Mostly, this is quite good, with only 3 outliers in the general case and one in
"Residues preceding proline". Hovering over any point in the plot will pop up a
label with the residue ID, and clicking on it will centre the view on that
residue. Find the outlier at Thr A84 (top right of the plot) and click on it.

.. figure:: images/3io0_Thr84.jpg

Hmm. Not only is this a Ramachandran outlier, but the red filled-in "cup"
between adjacent alpha carbons highlights that it's also a non-proline *cis*
peptide bond. These are exceedingly rare in real life, and only occur in
tightly-constrained, well-supported environments (which this isn't). It will
need to be fixed.

... But first, let's go through the other validation tools. Hide the
Ramachandran plot using the button at top right, and show the "Peptide bond
geometry" widget. This provides a list with all questionable (*cis* or twisted)
peptide bonds in the model. In this case we have seven non-proline *cis* bonds,
two severely twisted peptide bonds, and one *cis* proline (unlike the
non-proline case these are actually reasonably common at about 5% of proline
residues). Click through the list and have a look - while the *cis*-pro looks
happy, the rest appear unsupportable and will have to go.

Next, close the peptide bond widget and open the "Rotamers" widget. Click on
the first outlier in the list (Leu289).

.. figure:: images/3io0_Leu289.jpg

Just like with the backbone Ramachandran conformation, sidechain rotamer
outliers are marked up in real time with the exclamation mark/spiral motif you
see here. The size and colour of the marker will change depending on severity:
a large, pink one like this denotes a severe outlier, whereas a marginal
sidechain will show as a smaller yellow version.

Finally, close the rotamer widget and open the "Clashes" one. Once a simulation
has started (provided that all atoms are mobile) the molecular dynamics
forcefield ensures that clashes become effectively impossible - but the
presence of severe clashes in the starting coordinates means that we **will**
have to run a simulation of the entire model at least once to resolve them
before using more localised simulations to correct the errors we found above.

Let's go ahead and do that:

`isolde sim start #1`__

__ cxcmd:isolde\ sim\ start\ #1

If you're on a machine with a powerful GPU, this should only take a few seconds,
and a simulation running the entire model will likely be fast enough to work
with, without the need to stop and restart with smaller selections. On slower
hardware, just wait until you see the atoms start moving regularly - this
indicates that the minimisation is complete and clashes are resolved. Once you
hit the green stop button, it will be safe to start a simulation from any
smaller selection. The remainder of the tutorial will assume you're taking the
latter approach - if you are lucky enough to be working with a fast GPU, just
omit the stop/start steps.

Once your model is moving, take another look at the Clashes table - it should
now be reassuringly empty. Hide it, and go ahead and click the green stop
button. Now, let's take care of those *cis* bonds - since these are the most
"unusual" conformations in the model, they seem like a sensible place to start.
If you re-open the "Peptide bond geometry" widget and scroll to the bottom,
you'll see that the twisted peptides have taken care of themselves during the
energy minimisation - ISOLDE assumes that any peptide bond twisted more than 30
degrees from *cis* should probably actually be *trans*, and restrains it as
such.

Click on the top entry in the table, and switch to the "Rebuild" tab. You
should see that all the buttons in the top panel have now lighted up:

.. figure:: images/Thr84_remodel_options.png

The top two buttons provide the two possible types of peptide bond flip: a
simple 180 degree rotation (left), or a flip from *cis* to *trans* or *vice
versa* (right). We're currently interested in the latter.

**(NOTE: the buttons in this panel, like the others in the Rebuild tab, will
only be enabled when it is possible to use them. These two buttons, for example,
are only available when your selected atoms are from a single amino acid
residue with a peptide bond N-terminal to it. Similarly, the rotamer selection
buttons directly below them are only available when you have selected a
non-truncated amino acid residue with at least one sidechain chi angle.)**

Since we're not currently running a simulation, clicking either of these buttons
will automatically start a small localised simulation to perform the flip. Make
sure the text beside "Selected residue:" reads "A 84 THR" (if not, go back to
the Peptide bond geometry widget and click in the table again), then go ahead
and click the *cis<-->trans* button. After a few seconds initialisation time,
you should see something like this:

.. figure:: images/Thr84_pep_corrected.jpg

Voila! Not only is the *cis* peptide bond gone, but we've also corrected the
Ramachandran outlier that was here into the bargain. On the other hand, now the
sidechain is showing up as a severe rotamer outlier. Not to worry.

If you've deselected the residue, **ctrl-click** on an atom to select it again.
Now, the second-top row of buttons are the rotamer adjustment tools. The left
and right arrows scroll through previews of different rotamers (ordered by
statistical likelihood), while the following three buttons choose what to do
with the current preview: discard it, simply set it as the new coordinates, or
set it as (a) target angle(s) for torsion restraints to guide the simulation
smoothly there. If you add restraints, the fourth button will clear them. In
this case, the first rotamer in the list (accessed by clicking the right arrow)
is clearly the one we want - click that, then click the button to set the
coordinates to match the preview. Once you're happy with the results, click
the green "STOP" button.

Switch back to the "Validate" tab, and click the "Update list" button under the
peptide bond geometry table. Click the entry at the top of the table. If you
look around a little here, you'll see that there are actually *six* non-proline
*cis* bonds here, all in the same loop.

.. figure:: images/multi_cis_loop.jpg

Let's make a simulation selection that will cover all of them in one go. There
are two ways to go about this:

1. **ctrl-shift-click** on any atom to add it to the selection. Keep in mind
   that for each protein atom selected, when starting a simulation ISOLDE will
   automatically expand the selection to encompass its whole residue plus three
   residues forward and back along the chain.

2. Simply take advantage of ChimeraX's selection promotion framework. If you
   have individual atoms selected, pressing the up arrow once will expand the
   selection to whole residues, and up again will further expand it to whole
   secondary structure elements (helices, strands or unstructured stretches).
   Pressing the down arrow will reverse the process. In this case, after
   clicking on the entry in the Peptide bond geometry table the offending
   residue is selected. Make sure the focus is on the main GUI window by left-
   clicking anywhere on the molecule, then press up once to select the whole
   loop with all six *cis* bonds.

Once your selection is to your satisfaction, go ahead and click play to correct
them. You may choose to do this while everything is moving - or, if you prefer,
you can queue up all the flips with the simulation paused. Just **ctrl-click**
on a C-alpha C-terminal to a *cis* bond, then **ctrl-shift-click** on the
C-alphas for the other *cis* bonds you see. A handy way to remember which
peptide bond will flip is to keep in mind that if you've selected an alpha
carbon, the flip will be applied to the nitrogen directly bonded to it.
Now, type the command `isolde cisflip sel`__ in the command line. You will
not see a visual indicator that a flip of this type is pending, but you'll see
the results once you resume the simulation by pressing the play button. If you
accidentally flip one you didn't mean to, it's typically no big deal - just flip
it right back.

__ cxcmd:isolde\ cisflip\ sel


Once these flips are done, You will also most likely find that you need to
adjust the sidechains of Thr150 and Lys151 (although in some cases these correct
themselves once the *cis* bonds are fixed). Also keep in mind that the peptide
oxygen of Gly153 should end up pointing at the Arg148 sidechain - if this isn't
the case in your simulation, try using the "flip peptide bond" button (remember,
it's to the left of the *cis<-->trans* button).

**(NOTE: this form of peptide bond flip can also be applied using the command
"isolde pepflip sel". Both "isolde pepflip" and "isolde cisflip" should be used
with caution: while flipping every peptide bond in your model can certainly be
entertaining, it's not so much fun if you do so by accident.)**

Once you're done, the loop should look something like this:

.. figure:: images/multi_cis_loop_fixed.jpg

The remaining Ramachandran outlier (pink C-alpha) at Asp149 appears to be a
great example of the truism that not all Ramachandran outliers are wrong. In
this case, the conformation not only fits well to the density, but makes perfect
physical sense: the backbone twist is stabilised by charge interactions between
the carbonyl oxygens and two nearby arginine residues (Arg148 in the same chain
and Arg148 in the (y,z,x) symmetry copy), while the Asp149 sidechain is
salt-bridged to Arg120 and further stabilised by two hydrogen bonds. This is
what a Ramachandran outlier *should* look like: lots of energetically-favourable
interactions overcoming the energy penalty associated with the unusual backbone
twist.

You may have noticed some positive (green) difference density in the channel to
the right of this loop (when oriented like the image above), as well. This seems
quite persistent throughout refinement, and may represent some low-occupancy
ligand carried through during purification. At this resolution, though, it seems
impossible to identify without further experimental knowledge.

Anyway, once you stop the simulation and switch back to the "Validate" tab,
you'll find after clicking "Update list" that this has taken care of the last
of the non-proline *cis* peptides (and inspection of the remaining *cis* Pro
shows that it is indisputably real). So, we can hide the "Peptide bond geometry"
widget and move on.

From here, you could move on to working through the list of rotamer outliers in
a similar manner, but it remains a truism that in experimental model building,
human eyes should see each residue in context with its density at least once.
After all, it is entirely possible (and, in fact, quite common) for a model to
find an incorrect conformation that is *not* an outlier by any standard metric
but is nevertheless wrong.

The buttons at the far bottom right of the ISOLDE panel can help make this easy:

.. figure:: ../../../tools/gui/images/map_masking_toolbar.png

**(NOTE: these buttons are only available outside of simulations (that is, after
you've pressed a stop - not pause - button). During simulations the display is
set so that what you see matches exactly to what is simulated - while nothing
prevents you changing this via the command line, I'd strongly recommend against
it.)**

In particular, the left and right arrow buttons allow you to "step through" the
model in overlapping bite-sized chunks (specifically, two secondary structure
elements and their flanking loops), masking the map to the selection for easier
interpretation.  Clicking the right arrow once and zooming out should give you
a view like this:

.. figure:: images/structure_stepper.jpg

In my case, within the masked selection I see two rotamer outliers (Thr89 and
Val90) and one marginal rotamer (Thr79 - probably not so meaningful given its
high solvent exposure). However, at the first turn (between Gly87 and Asp88) I
also find this:

.. figure:: images/pep_flip_needed.jpg

This is a perfect example of why it's important to visually inspect the model.
Here we have no sidechain nor Ramachandran outliers, yet the telltale paired
red and green difference density blobs associated with the peptide oxygen and
nitrogen say that this conformation is unquestionably wrong. Perhaps more
concerning, in the original coordinates before energy minimisation in ISOLDE
there was not even any sign of trouble in the difference density, and only a
practised eye would notice the chemically-infeasible juxtaposition of the
peptide hydrogen with those of Arg120 and Val121:

.. figure:: images/pep_flip_needed_original_coords.jpg

This highlights what I consider to be one of the strengths of the approach
underlying ISOLDE: since the molecular dynamics forcefield is designed to
closely replicate the physical reality of the molecule, it is much less likely
to "fool you" in ways like this. Upon starting a simulation, the electrostatic
repulsion between these positively charged groups pushes the incorrectly
modelled peptide bond out of density, creating a clear red flag indicating
something is wrong. In this case, the remedy is clear: a simple peptide flip
creates two nice hydrogen bonds with the adjacent backbone (visualised here
using ChimeraX's H-bonds tool:

`hbonds #1`__

__ cxcmd:hbonds\ #1

.. figure:: images/pep_flipped.jpg

In addition, it places the backbone back in density (no more difference blobs),
and as a bonus converts the slightly-marginal Ramachandran case at Asp88 into a
favoured one. Don't forget to hide the H-bonds display:

`~hbonds`__

__ cxcmd:~hbonds


This covers most of the "workhorse" tools you're likely to use for day-to-day
work with models in medium-low resolution density with ISOLDE. Continue on
through the rest of the model at your leisure. After some practice, I've found
that I can improve this model to a MolProbity score below 1.0 in under 20
minutes (you can check for yourself using the MolProbity server_ or, if you have
it installed, *phenix.molprobity*). If you want to refine and/or check your work
using another package (recommended for any serious work: for a start, ISOLDE
currently does not attempt to refine atomic B-factors), remember you can save
the model and structure factors at any time using the following commands.

.. _server: http://molprobity.biochem.duke.edu/

`save ~/Desktop/isolde_tutorial.pdb #1`__

__ cxcmd:save\ ~/Desktop/isolde_tutorial.pdb\ #1

`save ~/Desktop/isolde_tutorial.mtz #1`__

__ cxcmd:save\ ~/Desktop/isolde_tutorial.mtz\ #1
