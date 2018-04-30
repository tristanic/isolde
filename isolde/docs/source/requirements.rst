System Requirements
===================

.. contents::
    :local:

Operating System
----------------

ISOLDE is currently distributed for Linux and MacOS operating systems, and has
essentially the same requirements as `ChimeraX`_. The most important requirement
in Linux is that the GCC runtime library must be version 4.9 and later. This is
true of Fedora 21 and above (although I have only personally tested version 24),
Ubuntu 16.04+, Debian 8. Additionally, ChimeraX and ISOLDE run happily in RedHat
7 environments (e.g. CentOS 7) if you provide the necessary libraries as
described `here <http://www.marlovitslab.org/chimerax/>`_. (*NOTE: rather than
copying the library files to the ChimeraX lib directory as described, I prefer
to place them in a separate directory and use a Bash script like the one below
to launch ChimeraX. This makes upgrading ChimeraX much easier.*)

.. code-block:: bash

  #!/bin/bash

  # Replace "/path/to" with the relevant paths on your machine as necessary

  CHIMERA_HOME=/path/to/chimerax

  export LD_LIBRARY_PATH=/path/to/gcc-4.9/lib64:/path/to/cuda/lib64:$LD_LIBRARY_PATH

  $CHIMERA_HOME/bin/ChimeraX $*


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
there are many similar laptops available) equipped with an NVidia GTX1070 GPU,
which allows me to work on >5,000-residue models and achieve interactive speeds
with simulations of up to 1,000 residues at a time. On the other end of the
scale, my 13" MacBook Air is a little trooper, achieving usable performance for
a few dozen residues at a time (enough to tweak a problematic loop or fix
local rotamer/backbone errors).





.. _ChimeraX: http://preview.cgl.ucsf.edu/chimerax/download.html
.. _Asus ROG Strix GL502VS: https://www.asus.com/us/ROG-Republic-Of-Gamers/ROG-GL502VS/
