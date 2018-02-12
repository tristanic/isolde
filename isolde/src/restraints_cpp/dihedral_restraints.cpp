
#include "dihedral_restraints.h"
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::Proper_Dihedral_Restraint_Mgr>;
template class pyinstance::PythonInstance<isolde::Proper_Dihedral_Restraint>;

namespace isolde
{

template <class DType, class RType>
isolde::Change_Tracker* Dihedral_Restraint_Base<DType, RType>::change_tracker() const
{
    return _mgr->change_tracker();
}

template <class DType, class RType>
RType* Dihedral_Restraint_Mgr_Base<DType, RType>::new_restraint(DType *d)
{
    auto it = _dihedral_to_restraint.find(d);
    if (it != _dihedral_to_restraint.end())
    {
        throw std::logic_error(error_duplicate());
        return nullptr;
    }
    return _new_restraint(d);
}

template <class DType, class RType>
RType* Dihedral_Restraint_Mgr_Base<DType, RType>::_new_restraint(DType *d)
{
    RType *r = new RType(d, this);
    _dihedral_to_restraint[d] = r;
    return r;
}

template <class DType, class RType>
RType* Dihedral_Restraint_Mgr_Base<DType, RType>::get_restraint(DType *d, bool create)
{
    auto it = _dihedral_to_restraint.find(d);
    if (it != _dihedral_to_restraint.end())
        return it->second;
    if (create)
        return _new_restraint(d);
    throw std::logic_error(error_no_restraint());
    return nullptr;
}


} // namespace isolde
