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

pkg_dir = "build/src/chimerax/clipper"  # DO NOT CHANGE

description = """
ChimeraX interface to Clipper libraries (http://www.ysbl.york.ac.uk/~cowtan/clipper/doc/)
for working with crystallographic data.
"""

pkg_data_list = ["/*.data", "*.so*", "*.dylib", "*.dll"]

for dir in "lib lib/clipper_python_core ".split():
   pkg_data_list.append(dir+'/*')

setup(
    name="ChimeraX-Clipper",
    version="0.1.0",  # PEP 440, should match Development Status below
    description="Interface to Kevin Cowtan's Clipper libraries",  # one line synopsis
    long_description=description,  # see above
    author="Tristan Croll",
    author_email="tic20@cam.ac.uk",
    url="",
    python_requires=">= 3.5",
    package_dir={
        "chimerax.clipper": pkg_dir,    # directory package's source files are in
    },
    packages=[
        "chimerax.clipper",
    ],
    package_data={
		"chimerax.clipper": pkg_data_list
	},
    install_requires=[
        "ChimeraX-Core >= 0.1",
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
        "License :: LGPL",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows :: Windows 10",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",

        "Development Status :: 2 - Pre-Alpha",
        "ChimeraX :: Bundle :: General :: 1,1 :: chimerax.clipper :: :: ",
        #"ChimeraX :: Tool :: Clipper :: General :: Clipper",
    ],
)
