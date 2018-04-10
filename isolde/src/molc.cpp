/*
 * Based upon ChimeraX molc.cpp. Interface between Python and C++
 * implementations of objects handling ChimeraX atomic data.
 */


#include <Python.h>     // Use PyUnicode_FromString

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>      // use PyArray_*(), NPY_*

#include <atomstruct/Atom.h>
#include <atomstruct/AtomicStructure.h>
#include <atomstruct/Bond.h>
#include <atomstruct/Chain.h>
#include <atomstruct/ChangeTracker.h>
#include <atomstruct/CoordSet.h>
#include <atomstruct/connect.h>
#include <atomstruct/destruct.h>     // Use DestructionObserver
#include <atomstruct/MolResId.h>
#include <atomstruct/PBGroup.h>
#include <atomstruct/polymer.h>
#include <atomstruct/Pseudobond.h>
#include <atomstruct/PBGroup.h>
#include <atomstruct/Residue.h>
#include <atomstruct/RibbonXSection.h>
#include <atomstruct/Ring.h>
#include <atomstruct/seq_assoc.h>
#include <atomstruct/Sequence.h>
#include <arrays/pythonarray.h>           // Use python_voidp_array()
#include <pysupport/convert.h>     // Use cset_of_chars_to_pyset

//#include "atomic_cpp/dihedral.h"
//#include "atomic_cpp/dihedral_mgr.h"
//#include "geometry/geometry.h"
//#include "interpolation/nd_interp.h"
//#include "validation/rama.h"
//#include "validation/rota.h"
//#include "restraints_cpp/changetracker.h"
//#include "restraints_cpp/position_restraints.h"
//#include "restraints_cpp/distance_restraints.h"
//#include "restraints_cpp/dihedral_restraints.h"
//#include "restraints_cpp/rotamer_restraints.h"
//#include "restraints_cpp/mdff.h"


#include <functional>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <stdint.h>
#include <string>
#include <vector>
#include <cmath>

#include "molc.h"

#include "atomic_cpp/dihedral_ext.h"
#include "atomic_cpp/dihedral_mgr_ext.h"

#include "validation/rama_ext.h"
#include "validation/rota_ext.h"

#include "restraints_cpp/changetracker_ext.h"
#include "restraints_cpp/mdff_ext.h"
#include "restraints_cpp/position_restraints_ext.h"
#include "restraints_cpp/distance_restraints_ext.h"
#include "restraints_cpp/dihedral_restraints_ext.h"
#include "restraints_cpp/rotamer_restraints_ext.h"
