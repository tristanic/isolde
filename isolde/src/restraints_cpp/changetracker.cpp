/**
 * @Author: Tristan Croll <tic20>
 * @Date:   26-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 28-Mar-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#define PYINSTANCE_EXPORT

#include "changetracker.h"
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::Change_Tracker>;

namespace isolde
{

Change_Tracker::Change_Tracker() {
    _reason_strings[REASON_RESTRAINT_CREATED] = std::string("created");
    _reason_strings[REASON_TARGET_CHANGED] = std::string("target changed");
    _reason_strings[REASON_SPRING_CONSTANT_CHANGED] = std::string("spring constant changed");
    _reason_strings[REASON_ADAPTIVE_C_CHANGED] = std::string("adaptive restraint constant changed");
    _reason_strings[REASON_DISPLAY_CHANGED] = std::string("display changed");
    _reason_strings[REASON_ENABLED_CHANGED] = std::string("enabled/disabled");
    _reason_strings[REASON_CUTOFF_CHANGED] = std::string("cutoff changed");
}

} //namespace isolde
