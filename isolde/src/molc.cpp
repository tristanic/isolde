/**
 * @Author: Tristan Croll <tic20>
 * @Date:   24-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 28-Mar-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



/*
 * Based upon ChimeraX molc.cpp. Interface between Python and C++
 * implementations of objects handling ChimeraX atomic data.
 */


#include <Python.h>     // Use PyUnicode_FromString

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>      // use PyArray_*(), NPY_*

#include <arrays/pythonarray.h>           // Use python_voidp_array()
#include <pysupport/convert.h>     // Use cset_of_chars_to_pyset

#include "molc.h"

#include "atomic_cpp/dihedral_ext.h"
#include "atomic_cpp/dihedral_mgr_ext.h"
#include "atomic_cpp/chiral_ext.h"
#include "atomic_cpp/chiral_mgr_ext.h"
#include "atomic_cpp/util.h"
#include "atomic_cpp/select.h"

#include "validation/rama_ext.h"
#include "validation/rota_ext.h"

#include "restraints_cpp/changetracker_ext.h"
#include "restraints_cpp/mdff_ext.h"
#include "restraints_cpp/position_restraints_ext.h"
#include "restraints_cpp/distance_restraints_ext.h"
#include "restraints_cpp/adaptive_distance_restraints_ext.h"
#include "restraints_cpp/dihedral_restraints_ext.h"
#include "restraints_cpp/rotamer_restraints_ext.h"
