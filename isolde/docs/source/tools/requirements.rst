.. _system_requirements:

System Requirements
===================

.. contents::
    :local:

Operating System
----------------

ISOLDE is distributed for Linux (CentOS 8 and later), MacOS and Windows 10
operating systems, and has essentially the same requirements as `ChimeraX`_.

Hardware
--------

The most important factor in getting the best out of ISOLDE is your GPU. ISOLDE
(or, to be more precise, the OpenMM library on which ISOLDE depends) relies on
GPU acceleration (via OpenCL or CUDA) to deliver sufficient molecular dynamics
performance for interactive simulations. As a minimum you will need an OpenCL
1.2 compatible GPU (don't worry - this describes almost every GPU manufactured
in at least the past 10 years), with the manufacturer's drivers installed.

ISOLDE is designed to provide optimum performance on consumer-level "gaming"
hardware. Currently, by far the best performance is seen on Nvidia GPUs.
Ideally, it should be run on a system with an Nvidia GTX 1060 or better. In the
later-generation RTX cards, the 2050 and 3050 should also give very good
performance, but their smaller on-board RAM may prove limiting when working with
very large models. On my current Windows test system (a laptop with Intel Core
i9-9900k CPU and a RTX 2080 GPU), a 1000 amino acid simulation runs with about 9
coordinate updates per second (at 50 simulation timesteps per coordinate update)
while maintaining >50 fps graphics. For non-interactive "settling" of big
models, ISOLDE successfully simulates 4v9o on a GTX 1070 GPU (8GB onboard RAM):
four (!) ribosomes in the asymmetric unit of a 2.9Å crystal - 991,533 atoms
including hydrogens. On the other end of the scale, my 13" MacBook Air is a
little trooper, achieving usable performance for a few dozen residues at a time
(enough to tweak a problematic loop or fix local rotamer/backbone errors). For
lower-end hardware like this, ISOLDE also provides a slightly lower-fidelity
"quick" simulation mode which provides significantly improved speed at the
expense of somehwat less realistic treatment of electrostatic interactions.





.. _ChimeraX: http://preview.cgl.ucsf.edu/chimerax/download.html
