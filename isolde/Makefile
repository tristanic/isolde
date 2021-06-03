# @Author: Tristan Croll <tic20>
# @Date:   28-Jan-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 02-Jan-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



# The "make" targets are:
# 	wheel: build a Python wheel in "dist" directory.
# 	app-install: build wheel (if needed) and install in ChimeraX.
# 	test: run ChimeraX
# 	debug: run ChimeraX with debugging flag set
# 	clean: remove files used in building wheel
# 	distclean: remove files used in building wheel and license file

# These parameters may be changed as needed.

# ChimeraX bundle names must start with "ChimeraX_"
# to avoid clashes with package names in pypi.python.org.
# When uploaded to the ChimeraX toolshed, the bundle
# will be displayed without the ChimeraX- prefix.
# BUNDLE_NAME = ChimeraX-ISOLDE
# BUNDLE_VERSION = 1.0b3
# ChimeraX bundles should only include packages
# that install as chimerax.package_name.
# General Python packages should be uploaded to
# pypi.python.org rather than the ChimeraX toolshed.
# PKG_NAME = chimerax.isolde

# Define where ChimeraX is installed.
OS = $(patsubst CYGWIN_NT%,CYGWIN_NT,$(shell uname -s))
# CHIMERAX_APP is the ChimeraX install folder
ifeq ($(OS),CYGWIN_NT)
ifndef RELEASE
# Windows
CHIMERAX_APP = "/c/Program Files/ChimeraX_Daily"
else
CHIMERAX_APP = "/c/Program Files/ChimeraX"
endif
endif

ifeq ($(OS),Darwin)
# Mac
ifndef RELEASE
CHIMERAX_APP = /Applications/ChimeraX_Daily.app
else
CHIMERAX_APP = /Applications/ChimeraX.app
endif
endif

ifeq ($(OS),Linux)
ifndef RELEASE
CHIMERAX_APP = /usr/libexec/UCSF-ChimeraX-daily
else
CHIMERAX_APP = /usr/libexec/UCSF-ChimeraX
endif
endif

# ==================================================================
# Theoretically, no changes are needed below this line

# Platform-dependent settings.  Should not need fixing.
# For Windows, we assume Cygwin is being used.
ifeq ($(OS),CYGWIN_NT)
CHIMERAX_EXE = $(CHIMERAX_APP)/bin/ChimeraX.exe
endif
ifeq ($(OS),Darwin)
CHIMERAX_EXE = $(CHIMERAX_APP)/Contents/bin/ChimeraX
export MACOSX_DEPLOYMENT_TARGET=10.13
endif
ifeq ($(OS),Linux)
CHIMERAX_EXE = $(CHIMERAX_APP)/bin/ChimeraX
endif

BUNDLE_BASE_NAME = $(subst ChimeraX-,,$(BUNDLE_NAME))
SOURCE = src
SRCS = $(SOURCE)/*.py $(SOURCE)/*.cpp

:DEFAULT_GOAL := wheel

#
# Actual make dependencies!
#



wheel $(WHEEL): bundle_info.xml $(SRCS)
	$(CHIMERAX_EXE) -m PyQt5.pyrcc_main -o src/resources/resources_rc.py src/resources/resources.qrc; \
	$(CHIMERAX_EXE) --nogui --safemode --cmd "devel build . ; exit"

install app-install: $(WHEEL)
	$(CHIMERAX_EXE) -m PyQt5.pyrcc_main -o src/resources/resources_rc.py src/resources/resources.qrc; \
	$(CHIMERAX_EXE) --nogui --safemode --cmd "devel install . ; exit"

uninstall app-uninstall:	$(WHEEL)
	$(CHIMERAX_EXE) --nogui --safemode --cmd "toolshed uninstall $(BUNDLE_BASE_NAME) ; exit"

docs:
	$(CHIMERAX_EXE) -m sphinx docs/source src/docs/user

test:
	$(CHIMERAX_EXE)

debug:
	$(CHIMERAX_EXE) --debug

clean:
	$(CHIMERAX_EXE) --nogui --safemode --cmd "devel clean . ; exit"

.PHONY: docs
