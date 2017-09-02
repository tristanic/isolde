# vim: set expandtab ts=4 sw=4:

# import distutils.core
# distutils.core.DEBUG = True
# TODO: remove distlib monkey patch when the wheel package
# implements PEP 426's pydist.json
from distlib import metadata
metadata.METADATA_FILENAME = "metadata.json"
from setuptools import setup
import os
import sys

pkg_dir = "build/src/chimerax/isolde"  # DO NOT CHANGE

description = """
ISOLDE: Interactive Structure Optimisation by Local Direct Exploration

ISOLDE is a next-generation environment for interactive building and
real-space refinement of atomic models into electron density maps. 
It applies interactive molecular dynamics to allow real-time, intuitive
performance of structural rearrangements from the small to the quite
drastic, while constantly maintaining physically reasonable interactions
with the surroundings.
"""

setup(
    name="ChimeraX-ISOLDE",
    version="0.9.13",  # PEP 440, should match Development Status below
    description="Interactive molecular dynamics based model building, manipulation and refinement",  # one line synopsis
    long_description=description,  # see above
    author="Tristan Croll",
    author_email="tic20@cam.ac.uk",
    url="http://www.ks.uiuc.edu/Research/mdff/",
    python_requires=">= 3.5",
    package_dir={
        "chimerax.isolde": pkg_dir,    # directory package's source files are in
    },
    packages=[
        "chimerax.isolde",
    ],
    package_data={
		"chimerax.isolde": ["molprobity_data/*.data", "amberff/*.xml", "amberff/amap/*.txt", "resources/*","*.so", "*.dylib", "*.dll"]
	},
    install_requires=[
        "ChimeraX-Core >= 0.1",
        #"ChimeraX-Clipper >= 0.1.0",
        # TODO: Should list numpy, OpenMM, simtk, etc.
    ],
    classifiers=[
        # From https://pypi.python.org/pypi?%3Aaction=list_classifiers
        # and our own ChimeraX-* classifiers.
        "Environment :: MacOS X :: Aqua",
        "Environment :: Win32 (MS Windows)",
        "Environment :: X11 Applications",
        "Framework :: ChimeraX",
        "Intended Audience :: Science/Research",
        "License :: Free for non-commercial use",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows :: Windows 10",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",

        "Development Status :: 2 - Pre-Alpha",
        "ChimeraX :: Bundle :: General :: 1,1 :: chimerax.isolde :: :: ",
        "ChimeraX :: Tool :: ISOLDE :: General :: Interactive Molecular Dynamics Flexible Fitting (iMDFF)",
        "ChimeraX :: Command :: fps :: General :: Display graphics FPS rate",
        "ChimeraX :: Command :: isolde :: General :: Command-line control of ISOLDE simulations",
    ],
)