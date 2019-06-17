System Requirements
===================

.. contents::
    :local:

Operating System
----------------

ISOLDE is distributed for Linux (CentOS 7 and later), MacOS and Windows 10
operating systems, and has essentially the same requirements as `ChimeraX`_.

Hardware
--------

The most important factor in getting the best out of ISOLDE is your GPU. ISOLDE
(or, to be more precise, the OpenMM library on which ISOLDE depends) relies on
GPU acceleration (via OpenCL or CUDA) to deliver sufficient molecular dynamics
performance for interactive simulations. As a minimum you will need an OpenCL
1.2 compatible GPU (don't worry - this describes almost every GPU manufactured
in at least the past 5 years), with the manufacturer's drivers installed.

ISOLDE is designed to provide optimum performance on consumer-level "gaming"
hardware. My regular test machine is an `Asus ROG Strix GL502VS`_ (not an ad -
there are many similar laptops available) equipped with an Nvidia GTX1070 GPU.
On this machine, "interactive" speeds (that is, better than about 10 fps) are
possible with over 1,000 residues mobile. For non-interactive "settling" of big
models, ISOLDE successfully simulates 4v9o: four (!) ribosomes in the
asymmetric unit of a 2.9Å crystal - 991,533 atoms including hydrogens. On the
other end of the scale, my 13" MacBook Air is a little trooper, achieving usable
performance for a few dozen residues at a time (enough to tweak a problematic
loop or fix local rotamer/backbone errors).





.. _ChimeraX: http://preview.cgl.ucsf.edu/chimerax/download.html
.. _Asus ROG Strix GL502VS: https://www.asus.com/us/ROG-Republic-Of-Gamers/ROG-GL502VS/
