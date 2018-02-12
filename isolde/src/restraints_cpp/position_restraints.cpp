
#define PYINSTANCE_EXPORT

#include "position_restraints.h"
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::Position_Restraint>;
template class pyinstance::PythonInstance<isolde::Position_Restraint_Mgr>;

namespace isolde
{

Change_Tracker* Position_Restraint::change_tracker() const
{
    return _mgr->change_tracker();
}

void Position_Restraint::target_vector(double *vector) const
{
    for (size_t i=0; i<3; ++i)
        *vector++ = _target[i]-_atom->coord()[i];
}

void Position_Restraint::set_k(double k)
{
    _spring_constant = k<0 ? 0.0 : ( k > MAX_SPRING_CONSTANT ? MAX_SPRING_CONSTANT : k);
}

double Position_Restraint::radius() const
{
    return _spring_constant/MAX_SPRING_CONSTANT
        * (RESTRAINT_MAX_RADIUS-RESTRAINT_MIN_RADIUS) + RESTRAINT_MIN_RADIUS;
}

void Position_Restraint::bond_cylinder_transform(float *rot44) const
{
    const Coord &xyz0 = atom()->coord();
    const Coord &xyz1 = get_target();
    double bvec[3];
    target_vector(bvec);
    auto &ac = atom()->coord();
    double d = geometry::l2_norm_3d(bvec);
    if (d == 0) {
        bvec[0]=0; bvec[1]=0; bvec[2]=0;
    } else {
        for (auto &b: bvec)
            b/=d;
    }
    double c = bvec[2], c1;
    if (c <= -1) c1 = 0;
    else c1 = 1.0/(1.0+c);
    double wx = -bvec[1], wy = bvec[0];
    double cx = c1*wx, cy = c1*wy;
    double r = radius();
    double h = d*2;
    *rot44++ = r*(cx*wx + c);
    *rot44++ = r*cy*wx;
    *rot44++ = -r*wy;
    *rot44++ = 0;

    *rot44++ = r*cx*wy;
    *rot44++ = r*(cy*wy + c);
    *rot44++ = r*wx;
    *rot44++ = 0;

    *rot44++ = h*wy;
    *rot44++ = -h*wx;
    *rot44++ = h*c;
    *rot44++ = 0;

    *rot44++ = (xyz0[0]+xyz1[0])/2;
    *rot44++ = (xyz0[1]+xyz1[1])/2;
    *rot44++ = (xyz0[2]+xyz1[2])/2;
    *rot44++ = 1;

}

Position_Restraint* Position_Restraint_Mgr::_new_restraint(Atom *atom, const Coord& target)
{
    if (atom->structure() != _atomic_model) {
        throw std::logic_error(error_different_mol());
    }
    if (atom->element().number() == 1) {
        throw std::logic_error(error_hydrogen());
    }
    Position_Restraint* restraint = new Position_Restraint(atom, target, this);
    _atom_to_restraint[atom] = restraint;
    return restraint;
}

Position_Restraint* Position_Restraint_Mgr::_new_restraint(Atom *atom)
{
    return _new_restraint(atom, atom->coord());
}

Position_Restraint* Position_Restraint_Mgr::get_restraint(Atom *atom, bool create)
{
    auto it = _atom_to_restraint.find(atom);
    if (it != _atom_to_restraint.end())
        return it->second;
    if (create)
        return _new_restraint(atom);
    return nullptr;
}

std::vector<Position_Restraint *> Position_Restraint_Mgr::visible_restraints() const
{
    std::vector<Position_Restraint *> visibles;
    for (auto &it: _atom_to_restraint)
    {
        auto r = it.second;
        if (r->visible())
            visibles.push_back(r);
    }
    return visibles;
}


void Position_Restraint_Mgr::delete_restraints(const std::set<Position_Restraint *>& to_delete)
{
    auto db = DestructionBatcher(this);
    _delete_restraints(to_delete);
}

void Position_Restraint_Mgr::_delete_restraints(const std::set<Position_Restraint *>& to_delete)
{
    for (auto r: to_delete) {
        _atom_to_restraint.erase(r->atom());
        delete r;
    }
}

void Position_Restraint_Mgr::destructors_done(const std::set<void *>& destroyed)
{
    auto db = DestructionBatcher(this);
    std::set<Position_Restraint *> to_delete;
    for (auto &it: _atom_to_restraint) {
        auto r = it.second;
        if (destroyed.find(static_cast<void *>(r->atom())) != destroyed.end()) {
            to_delete.insert(r);
        }
    }
    _delete_restraints(to_delete);
}

Position_Restraint_Mgr::~Position_Restraint_Mgr()
{
    auto du = DestructionUser(this);
    for (auto &it: _atom_to_restraint) {
        delete it.second;
    }
}

} //namespace isolde
