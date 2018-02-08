#ifndef ISOLDE_POSITION_RESTRAINTS
#define ISOLDE_POSITION_RESTRAINTS

#include <atomstruct/destruct.h>
#include <atomstruct/string_types.h>
#include <atomstruct/Atom.h>
#include <atomstruct/Coord.h>
#include <atomstruct/Pseudobond.h>
#include <atomstruct/Residue.h>
#include <pyinstance/PythonInstance.declare.h>

using namespace atomstruct;

namespace isolde
{
class Position_Restraint_Mgr;

class Position_Restraint: public pyinstance::PythonInstance<Position_Restraint>
{
public:
    Position_Restraint() {}
    ~Position_Restraint() { auto du=DestructionUser(this); }
    Position_Restraint(Atom* atom, const Coord& target)
        : _atom(atom)
    {
        set_target(target[0], target[1], target[2]);
    }

    void set_target(const Real &x, const Real &y, const Real &z)
    {
        _target[0]=x; _target[1]=y; _target[2]=z;
    }
    void set_target(Real *target)
    {
        for (size_t i=0; i<3; ++i)
            _target[i] = *(target++);
    }
    const Coord& get_target() const { return _target; }
    void get_target(double *target) const {
        for (size_t i=0; i<3; ++i)
            *target++ = _target[i];
    }
    void set_spring_constant(const double &k) { _spring_constant = k; }
    double spring_constant() const { return _spring_constant; }
    void set_enabled(bool flag) { _enabled = flag; }
    bool enabled() const { return _enabled; }
    bool visible() const { return _atom->visible() && _enabled; }
    void target_vector(double *vector) const;
    Atom *atom() const { return _atom; }

private:
    Atom* _atom;
    Coord _target;
    double _spring_constant = 0.0;
    bool _enabled = false;


}; //class Position_Restraint

class Position_Restraint_Mgr
    : public DestructionObserver, public pyinstance::PythonInstance<Position_Restraint_Mgr>
{
public:
    Position_Restraint_Mgr() {}
    ~Position_Restraint_Mgr();
    Position_Restraint_Mgr(Structure *atomic_model)
        : _atomic_model(atomic_model) {}

    Structure* structure() const { return _atomic_model; }
    Position_Restraint* get_restraint(Atom *atom, bool create);
    size_t num_restraints() const { return _atom_to_restraint.size(); }

    void delete_restraints(const std::set<Position_Restraint *>& to_delete);
    virtual void destructors_done(const std::set<void *>& destroyed);

private:
    Structure* _atomic_model;
    std::unordered_map<Atom*, Position_Restraint*> _atom_to_restraint;
    Position_Restraint* _new_restraint(Atom *atom);
    Position_Restraint* _new_restraint(Atom *atom, const Coord& target);
    const char* error_different_mol() const {
        return "This atom is in the wrong structure!";
    }
    const char* error_hydrogen() const {
        return "Restraints on hydrogen atoms are not allowed!";
    }
    void _delete_restraints(const std::set<Position_Restraint *>& to_delete);

}; //class Position_Restraint_Mgr

} //namespace isolde
#endif //ISOLDE_POSITION_RESTRAINTS
