![ISOLDE logo](https://github.com/tristanic/isolde/blob/master/logo/isolde_logo.jpg)

# ISOLDE
Interactive molecular dynamics based model building into low-resolution crystallographic and cryo-EM maps

[Home page](https://isolde.cimr.cam.ac.uk/)

> **⚠️ `rdkit` branch — work in progress.**
> This branch adds an RDKit-based chemistry layer (chirality-aware rebuilding,
> automatic chiral-restraint generation, and ligand registration/placement). It
> depends on the **`ChimeraX-ChemComp`** bundle, which is **not yet publicly
> available**, so it cannot currently be built or installed outside Altos Labs.
> The dependency is isolated behind a small seam and can be dropped (or
> `ChemComp` vendored into ISOLDE) once it is released. **For general use, build
> from the [`master`](https://github.com/tristanic/isolde/tree/master) branch.**

> **⚠️ `garnet-ff` branch — also work in progress.**
> On top of the `rdkit` layer above, this branch integrates the experimental
> **GARNET** graph-ML force field as an opt-in alternative to AMBER (selected with
> `isolde set forcefield garnet`; AMBER stays the default). It needs a few extra
> Python packages and the `garnet_core` source tree, which are not part of a plain
> ISOLDE install — see [Experimental: the GARNET force field](#experimental-the-garnet-force-field-development)
> below.

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

In Windows, a special executable is needed for the command-line switches to be correctly handled. The following command will clean, build and install ISOLDE:

`C:\Program Files\ChimeraX\bin\ChimeraX-console.exe" -m PyQt5.pyrcc_main -o src/resources/resources_rc.py src/resources/resources.qrc`

`C:\Program Files\ChimeraX\bin\ChimeraX-console.exe" --nogui --cmd "devel clean .; devel install .; exit"`

For convenience, these commands are also wrapped in a batch file, make_win.bat.

`make_win release clean app-install`

... will achieve the same result as the above.

## Experimental: the GARNET force field (development)

> Status: development-grade and opt-in. AMBER remains the default; nothing here
> runs unless you explicitly select the GARNET force field. There is deliberately
> no polished install pipeline yet.

On the `garnet-ff` branch, an ISOLDE simulation can be parameterised by the
[garnet-isolde](https://github.com/altos-labs/garnet-isolde) graph-ML force field
instead of by AMBER template matching. Two things are needed beyond a normal
ISOLDE build, both installed into **ChimeraX's own Python**:

1. **Python packages** — `torch` and `torch-geometric`. These are *not* ISOLDE
   dependencies (an AMBER-only install should not pay a ~1 GB torch download), and
   nothing in ISOLDE proper imports them. They are declared by the `garnet-isolde`
   distribution in step 2, so installing that pulls them in automatically. Two
   caveats:
   - `torch` must be **>= 2.13** (what the wheel pins, matching the interpreter the
     shipped checkpoints were trained under). Anything <= 2.9 has a NaN in the
     `atan2(0, 0)` gradient that GARNET's dihedral term can hit. Note that torch
     2.8 and 2.11+ also disagree on that backward, so *relaxed geometries* are
     torch-version dependent even when both work; forward energies are not.
   - On a machine with an NVIDIA GPU, install the matching **CUDA** wheel of
     `torch` manually **first** (from https://pytorch.org), before installing
     `garnet-isolde`, so pip doesn't settle on the CPU-only build.

2. **`garnet_core` and the trained weights** — from the `garnet-isolde`
   repository. It is **not on PyPI**, so it is installed by hand into ChimeraX's
   Python. Two ways, depending on whether you develop the force field or just use it.

   **If you were handed a wheel** (the usual case — it bundles the trained
   checkpoints, and is pure Python, so the one file works on Windows, Linux and
   macOS):

   Linux/macOS:
   ```
   /path/to/ChimeraX/bin/ChimeraX --nogui --cmd "pip install garnet-isolde@file:///path/to/garnet_isolde-<version>-py3-none-any.whl ; exit"
   ```
   Windows:
   ```
   "C:\Program Files\ChimeraX\bin\ChimeraX-console.exe" --nogui --cmd "pip install garnet-isolde@file:///C:/path/to/garnet_isolde-<version>-py3-none-any.whl ; exit"
   ```

   > Note the `garnet-isolde@file://...` form, and the *three* slashes. ChimeraX's
   > `pip` command validates its argument as a PEP 508 requirement and rejects a bare
   > `.whl` path with *"invalid requirement specified"*; a direct reference is a valid
   > requirement, so this is the way through. (Alternatively, bypass the ChimeraX
   > command entirely: `PYTHONUSERBASE=<ChimeraX user dir> <ChimeraX>/bin/python3.x -m
   > pip install --user <wheel>`.)

   To build that wheel from a checkout, run `pip wheel --no-deps -w dist .` at the
   repo root with any Python >= 3.11; the result lands in `dist/`.

   **If you develop `garnet_core` itself**, install the checkout *editable* instead,
   so your edits take effect without reinstalling:
   ```
   /path/to/ChimeraX/bin/ChimeraX --nogui --cmd "pip install -e /path/to/garnet-isolde ; exit"
   ```

   Either way the bundled checkpoints are found automatically: ISOLDE pins each
   variant by bare filename and `garnet_core.weights` resolves it, knowing both the
   installed layout (`garnet_core/trained_models/`) and a checkout's
   `garnetff/trained_models/`. To point the bare `garnet` alias at a one-off
   checkpoint, set the `ISOLDE_GARNET_CHECKPOINT` environment variable to its path.

   > The wheel is an internal build: the licence for this fork of GARNET is not yet
   > settled, so it carries the `Private :: Do Not Upload` marker and is not for
   > distribution outside the team. See its bundled `NOTICE`.

Once both are in place, start ISOLDE, set the experience level to **Developer**
(a force-field selector then appears in ISOLDE's *General* tab), or simply run
`isolde set forcefield <name>`. Switch back to AMBER with
`isolde set forcefield amber14`.

GARNET is offered as one entry per training run, named `garnet-{run}`, so different
incarnations can be selected and compared side by side in the one session — each pinned
to its own checkpoint, with its OpenMM functional form auto-detected from that checkpoint
(later rounds change the form, e.g. `garnet-r10b` adds a per-atom repulsive wall and a
short-range Coulomb guard that `garnet-r5d` lacks). Currently available:

| name | checkpoint | functional form |
|---|---|---|
| `garnet-r5d`  | `dtr_sf_r5d_ep1.pt`  | double-exponential vdW, scalar wall exponent |
| `garnet-r10b` | `dtr_sf_r10b_ep2.pt` | per-atom repulsive wall + short-range Coulomb guard |

The bare `garnet` name is kept as a backward-compatible alias for the default checkpoint.
Add a future round by dropping its `garnet-{run}` entry into `_GARNET_VARIANTS`
(`isolde/src/openmm/forcefields.py`) and a matching profile — no other code changes.

The planned chemistry-verification framework that will accompany this
(connecting every component to a verified CCD/SMILES source of truth, and
handling incomplete models) is described in
`isolde/src/openmm/garnet/CHEMISTRY_PROVENANCE_BRIEF.md`.

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

`C:\Program Files\ChimeraX\bin\ChimeraX-console.exe" -m sphinx docs/source src/docs/user`

or

`make_docs.bat release`

... then reinstall ISOLDE.

Remember, you can always find the most recent version of the documentation [here](https://isolde.cimr.cam.ac.uk/documentation/).
