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
}

} //namespace isolde
