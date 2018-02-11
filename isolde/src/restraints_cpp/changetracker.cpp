#define PYINSTANCE_EXPORT

#include "changetracker.h"
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::Change_Tracker>;

namespace isolde
{

Change_Tracker::Change_Tracker() {
    _reason_strings[REASON_DISTANCE_RESTRAINT_CREATED] = std::string("Distance restraint created");
    _reason_strings[REASON_DISTANCE_RESTRAINT_CHANGED] = std::string("Distance restraint changed");
    _reason_strings[REASON_DIHEDRAL_RESTRAINT_CREATED] = std::string("Dihedral restraint created");
    _reason_strings[REASON_DIHEDRAL_RESTRAINT_CHANGED] = std::string("Dihedral restraint changed");
    _reason_strings[REASON_POSITION_RESTRAINT_CREATED] = std::string("Position restraint created");
    _reason_strings[REASON_POSITION_RESTRAINT_CHANGED] = std::string("Position restraint changed");
}

} //namespace isolde
