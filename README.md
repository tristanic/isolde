![ISOLDE logo](https://github.com/tristanic/isolde/blob/master/logo/isolde_logo.jpg)

# ISOLDE
Interactive molecular dynamics based model building into low-resolution crystallographic and cryo-EM maps

[Home page](https://isolde.cimr.cam.ac.uk/)

## What is ISOLDE?

ISOLDE is a plugin to [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/), designed to ease the task of
macromolecular model building into low-to-medium resolution maps derived from crystallographic or electron
cryo-microscopy (cryoEM) experiments. Historically, this has been an extremely challenging task, since at
low resolutions the maps alone are insufficient to precisely place individual atoms. ISOLDE aims to reduce
this challenge in a number of ways:

* Rebuilding is accomplished via GPU-accelerated interactive molecular dynamics (using [OpenMM](http://openmm.org/)
  and the [AMBER molecular dynamics forcefield](https://ambermd.org/AmberModels.php)) to make the task feel as close
  as possible to what it might be like to work with a real physical molecule.
* Geometric validation of protein backbone and sidechain conformations is performed in real time, allowing you to see
  problem sites directly on the model as you work with it.
* Remodelling can be performed by directly tugging on atoms, or via the interactive addition and removal of position,
  torsion and/or distance restraints
* For crystallographic datasets, structure factors are constantly recalculated in the background as the model coordinates
  change - as the model improves, you see the map improve.

## What does it look like?

Like this:

![ISOLDE example image](https://github.com/tristanic/isolde/blob/master/isolde/docs/source/tutorials/intro/crystal_intro/images/3io0_Thr84.jpg)

For other examples and demonstration videos, see [the ISOLDE webpage](https://isolde.cimr.cam.ac.uk).

## How do I get it?

In most cases, you should not need to build ISOLDE from source for yourself. Regular (approximately fortnightly)
development builds are released for Linux, Mac and Windows on the ChimeraX Tool Shed, and can be installed
directly from within ChimeraX itself. In general, just download and install the latest daily build of ChimeraX
from [here](https://www.cgl.ucsf.edu/chimerax/download.html#daily), then go to Tools/More Tools... and follow
the links to ISOLDE.

## Compiling from source

**NOTE:** Some large files in this repository are stored using [Git-LFS](https://git-lfs.github.com/). To clone these to your own system you'll need to have the Git-LFS client installed.

ISOLDE uses ChimeraX's [bundle building pipeline](https://www.cgl.ucsf.edu/chimerax/docs/devel/writing_bundles.html), with
the majority of the build information defined in [bundle_info.xml](https://github.com/tristanic/isolde/blob/master/isolde/bundle_info.xml).
Dependencies outside of those already present in ChimeraX itself are kept to a minimum. You will need to have a compatible
version of ISOLDE's sister package [ChimeraX-Clipper](https://github.com/tristanic/chimerax-clipper) installed (whether built
from source or installed via the Tool Shed). Additionally, you will need to provide the OpenMM header files (these are not
currently distributed with ChimeraX). The paths to these are hard-coded in bundle_info.xml, so you'll need to change the
following lines to the correct path(s):

```xml
      <IncludeDir platform="mac">/Users/tic20/anaconda3/envs/openmm74/include</IncludeDir>
      <IncludeDir platform="linux">/home/tic20/anaconda3/envs/openmm74/include</IncludeDir>
      <IncludeDir platform="windows">C:\Users\tic20\Anaconda3\envs\openmm74\include</IncludeDir>
```

**IMPORTANT NOTE FOR LINUX USERS**: the version of GCC you use for building needs to be binary-compatible with the version
used to build **both** ChimeraX and OpenMM. The version of OpenMM distributed with ChimeraX is the official release, compiled
with GCC 4.8. The RedHat and Generic Linux builds of ChimeraX are both built with GCC 4.9, so your easiest path is to also
build ISOLDE using GCC 4.9 (if you're using a RedHat flavour, the most convenient way is to use devtoolset-3). The Ubuntu
builds of ChimeraX are built with the GCC versions shipped with the OS - these are *not* backward-compatible. If you really
wish to build with these, you will first need to build and install your own version of OpenMM into the ChimeraX environment
following the instructions [here](http://docs.openmm.org/latest/userguide/library.html#compiling-openmm-from-source-code).

For MacOS you will need to have XCode installed, while for Windows you will need Visual Studio 2015 or better.

Once the above conditions are met, you can go ahead and build in a Linux or MacOS environment as follows:

- change to the directory containing bundle_info.xml
- run the following commands:

`/path/to/ChimeraX/bin/ChimeraX -m PyQt5.pyrcc_main -o src/resources/resources_rc.py src/resources/resources.qrc`

`/path/to/ChimeraX/bin/ChimeraX --nogui --cmd "devel build . ; exit"`

- then install with:

`/path/to/ChimeraX/bin/ChimeraX --nogui --cmd "devel install . ; exit"`

To clean the compilation and start from scratch, use:

`/path/to/ChimeraX/bin/ChimeraX --nogui --cmd "devel clean . ; exit"`

For convenience, these are also wrapped in a simple Makefile, allowing the above to be achieved with the equivalent commands:

`make`, `make install` and `make clean`

respectively.

In Windows, the corresponding commands are slighty more complex - an extra layer of indirection is needed for the
command-line switches to be correctly handled. The following command will clean, build and install ISOLDE:

`C:\Program Files\ChimeraX\bin\python.exe" "c:\Program Files\ChimeraX\bin\Lib\site-packages\ChimeraX_main.py" -m PyQt5.pyrcc_main -o src/resources/resources_rc.py src/resources/resources.qrc`

`C:\Program Files\ChimeraX\bin\python.exe" "c:\Program Files\ChimeraX\bin\Lib\site-packages\ChimeraX_main.py" --nogui --cmd "devel clean .; devel install .; exit"`

## Building the documentation

ISOLDE's documentation tree is defined using [Sphinx](http://www.sphinx-doc.org/en/master/). While Sphinx itself is bundled
with ChimeraX, you will also need an installation of LaTeX on your system path.

Since Sphinx's source code documentation relies on introspection from within Python itself, you will need to have already built
and installed ISOLDE into ChimeraX before the documentation can be built. Once you've done that, change to the directory
containing bundle_info.xml, then do the following:

Linux/MacOS:

`/path/to/ChimeraX/bin/ChimeraX -m sphinx docs/source src/docs/user`

or

`make docs`

Windows:

`C:\Program Files\ChimeraX\bin\python.exe" "c:\Program Files\ChimeraX\bin\Lib\site-packages\ChimeraX_main.py" -m sphinx docs/source src/docs/user`

... then reinstall ISOLDE.

Remember, you can always find the most recent version of the documentation [here](https://isolde.cimr.cam.ac.uk/documentation/).
