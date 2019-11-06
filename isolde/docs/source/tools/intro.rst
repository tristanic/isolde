What is ISOLDE all about?
=========================

Simply put, the goal of ISOLDE is to facilitate the building of high-quality
atomic models into moderate-to-low resolution experimental maps, where the
experimental information alone is insufficient to precisely place individual
atoms. Historically this has been an exceedingly challenging task, for a
number of reasons:

    (a) our understanding of how real atoms and molecules behave was limited;
    (b) we didn't have the computational resources to realistically model them
        at useful speeds anyway; and
    (c) computer graphics were not up to the task of providing intelligible
        and information-rich interactive visualisations.

Today, none of the above remains true. Molecular dynamics forcefields such as
AMBER and CHARMM provide high-fidelity descriptions of the forces governing
most macromolecules and a growing population of small molecule ligands.
Molecular dynamics engines such as OpenMM leverage the massively-parallel
computing capabilities of modern graphics processing units (GPUs) to solve
Newton's laws of motion for these forces hundreds of times per second on
systems of a few thousand atoms. Meanwhile, ChimeraX provides a fast, flexible
API allowing the clear, high-speed and rich rendering of the ongoing simulation
necessary for it to make sense to human eyes.

Combining the above allows the model-building task to be re-imagined as a truly
interactive experience in which, rather than carefully adjusting individual
atoms and dihedral angles, the user instead helps to guide a "living", explicitly
physical model into the experimental map. In this way, many of the unlikely or
impossible atomic arrangements that used to plague model builders are prevented
from ever occurring in the first place, and many other errors simply fix
themselves.

Alongside the above, a core component of ISOLDE's design philosophy is the
need for real-time, continuous feedback. In ISOLDE the process of rotamer,
Ramachandran and peptide plane validation have been streamlined to bring their
combined time  down to 1-2 milliseconds for a 500-residue structure. This allows
these core validation metrics to be evaluated in real time every time the model
coordinates change, with visualisations showing you at a glance exactly which
residues are problematic at any moment in time. Expect to see further live
validation features appearing as ISOLDE grows.
