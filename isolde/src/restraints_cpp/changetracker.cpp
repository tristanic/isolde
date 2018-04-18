/**
 * @Author: Tristan Croll
 * @Date:   13-Feb-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   Tristan Croll
 * @Last modified time: 18-Apr-2018
 * @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
 * @Copyright: Copyright 2017-2018 Tristan Croll
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
    _reason_strings[REASON_DISPLAY_CHANGED] = std::string("display changed");
    _reason_strings[REASON_ENABLED_CHANGED] = std::string("enabled/disabled");
    _reason_strings[REASON_CUTOFF_CHANGED] = std::string("cutoff changed");
}

} //namespace isolde
