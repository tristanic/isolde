/**
 * @Author: Tristan Croll <tic20>
 * @Date:   26-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 02-Apr-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#define PYINSTANCE_EXPORT
#include "dihedral_restraints.h"
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::Chiral_Restraint_Mgr>;
template class pyinstance::PythonInstance<isolde::Chiral_Restraint>;

template class pyinstance::PythonInstance<isolde::Proper_Dihedral_Restraint_Mgr>;
template class pyinstance::PythonInstance<isolde::Proper_Dihedral_Restraint>;

namespace isolde
{

template <class DType, class RType>
void Dihedral_Restraint_Mgr_Base<DType, RType>::_delete_restraints(const std::set<RType *>& to_delete)
{
    for (auto r: to_delete) {
        auto d = r->get_dihedral();
        _dihedral_to_restraint.erase(d);
        delete r;
    }
}

template <class DType, class RType>
void Dihedral_Restraint_Mgr_Base<DType, RType>::delete_restraints(const std::set<RType *>& to_delete)
{
    auto db = DestructionBatcher(this);
    _delete_restraints(to_delete);
}

template <class DType, class RType>
void Dihedral_Restraint_Mgr_Base<DType, RType>::destructors_done(const std::set<void *>& destroyed)
{
    auto db = DestructionBatcher(this);
    std::set<RType *> to_delete;
    for (auto &it: _dihedral_to_restraint) {
        auto d = it.first;
        if (destroyed.find(d) != destroyed.end())
            to_delete.insert(it.second);
    }
    _delete_restraints(to_delete);
}

template <class DType, class RType>
Dihedral_Restraint_Mgr_Base<DType, RType>::~Dihedral_Restraint_Mgr_Base()
{
    auto du = DestructionUser(this);
    for (auto &it: _dihedral_to_restraint)
        delete it.second;
    _dihedral_to_restraint.clear();
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
    RType *r = new RType(d, static_cast<Dihedral_Restraint_Change_Mgr *>(this));
    _dihedral_to_restraint[d] = r;
    track_created(static_cast<const void *>(r));
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
    // throw std::logic_error(error_no_restraint());
    return nullptr;
}

template <class DType, class RType>
std::vector<RType *> Dihedral_Restraint_Mgr_Base<DType, RType>::visible_restraints() const
{
    std::vector<RType *> visibles;
    for (auto &it: _dihedral_to_restraint)
    {
        auto r = it.second;
        if (r->visible())
            visibles.push_back(r);
    }
    return visibles;
}

/***************************************************
 *
 * Specialisations
 *
 ***************************************************/

Chiral_Restraint::Chiral_Restraint(
    Chiral_Center *chiral, Dihedral_Restraint_Change_Mgr *mgr)
    : Dihedral_Restraint_Base<Chiral_Center>(chiral, mgr)
    {
        // Chiral_Restraint instances are enabled by default.
        _enabled = true;
        _spring_constant = DEFAULT_CHIRAL_RESTRAINT_SPRING_CONSTANT;
        _cutoff = DEFAULT_CHIRAL_RESTRAINT_CUTOFF;
    }


 Proper_Dihedral_Restraint::Proper_Dihedral_Restraint(
     Proper_Dihedral *dihedral, Dihedral_Restraint_Change_Mgr *mgr)
     : Dihedral_Restraint_Base<Proper_Dihedral>(dihedral, mgr)
     {}


//! Provide the transforms for the Proper_Dihedral_Restraint annotation
/*! The Proper_Dihedral_Restraint annotation requires *two* transforms:
 *  one for the rotation of the "ring" relative to the "post", and one to
 *  scale and move the result to the position of the axial bond.
 *  The transform matrices are in 3x4 format where the left-hand-side
 *  3x3 is the rotation and the right-hand-side is the translation.
 */
void Proper_Dihedral_Restraint::get_annotation_transform(float *tf)
{
    // First the rotational component (rotating the ring relative to the post)
    // Sets the direction the arrow points
    float *tf1 = tf;
    float *tf2 = tf+16;
    bool flip = offset() < 0;
    float width = get_spring_constant()/MAX_RADIAL_SPRING_CONSTANT *
        (DIHEDRAL_RESTRAINT_MAX_WIDTH-DIHEDRAL_RESTRAINT_MIN_WIDTH) + DIHEDRAL_RESTRAINT_MIN_WIDTH;

    geometry::rotation_gl<float>(Z_AXIS, offset(), tf1);
    if(flip)
    // Flipping on x axis - negate y and z components
    {
        geometry::flip_on_x_gl(tf1);
    }
    auto b = get_dihedral()->axial_bond();
    const auto &c0 = b->atoms()[0]->coord();
    const auto &c1 = b->atoms()[1]->coord();
    // float xyz0[3], xyz1[3];
    // for (size_t i=0; i<3; ++i)
    // {
    //     xyz0[i]=c0[i];
    //     xyz1[i]=c1[i];
    // }

    // get the transformation mapping to the bond position
    geometry::bond_cylinder_transform_gl<Coord, float>(c0, c1, 1.0, width, tf2);
    // apply the bond transformation to the rotational component
    float temp[16];
    geometry::multiply_transforms_gl<float>(tf1, tf2, temp);
    for (size_t i = 0; i<16; ++i) {
        tf1[i] = temp[i];
    }
} //get_annontation_transform

void Proper_Dihedral_Restraint::get_annotation_color(uint8_t *color)
{
    colors::color thecolor;
    colormap()->interpolate(std::abs(offset()), _cutoffs, thecolor);
    for (size_t j=0; j<4; ++j) {
        *(color++) = (uint8_t)(thecolor[j]*255.0);
    }
} //get_annotation_color


template class Dihedral_Restraint_Base<Chiral_Center>;
template class Dihedral_Restraint_Mgr_Base<Chiral_Center, Chiral_Restraint>;

template class Dihedral_Restraint_Base<Proper_Dihedral>;
template class Dihedral_Restraint_Mgr_Base<Proper_Dihedral, Proper_Dihedral_Restraint>;

} // namespace isolde
