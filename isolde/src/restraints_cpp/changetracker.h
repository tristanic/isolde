#ifndef ISOLDE_CHANGETRACKER
#define ISOLDE_CHANGETRACKER

#include <vector>
#include <string>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <set>

#include <pyinstance/PythonInstance.declare.h>

namespace isolde
{

class Change_Tracker: public pyinstance::PythonInstance<Change_Tracker>
{
public:
    typedef std::unordered_map<const void*, std::set<const void *>> mgr_to_managed;
    typedef std::unordered_map<std::type_index, mgr_to_changeds> ptr_type_to_set;

    enum Reasons{
        REASON_DISTANCE_RESTRAINT_CREATED,
        REASON_DISTANCE_RESTRAINT_CHANGED,
        REASON_DIHEDRAL_RESTRAINT_CREATED,
        REASON_DIHEDRAL_RESTRAINT_CHANGED,
        REASON_POSITION_RESTRAINT_CREATED,
        REASON_POSITION_RESTRAINT_CHANGED
    };

    Change_Tracker();
    ~Change_Tracker() {}

    template <class Mgr, class C>
    void add_created(Mgr* mgr, C* ptr)
    {
        auto &m_created = _created_map[std::type_index(type(mgr))];
        m_created[mgr].insert(ptr);
    }

    template <class Mgr, class C>
    void add_created_batch(Mgr* mgr, std::vector<C*> ptrs)
    {
        auto &m_created = _created_map[std::type_index(type(mgr))];
        auto &created = m_created[mgr];
        for (const auto &p: ptrs)
            created.insert(p);
    }

    template <class Mgr, class C>
    void add_modified(Mgr *mgr, C* ptr)
    {
        auto &m_changes = _changed_map[std::type_index(type(mgr))];
        m_changes[mgr].insert(ptr);
    }

    template <class Mgr, class C>
    void add_modified_batch(Mgr *mgr, std::vector<C*> ptrs)
    {
        auto &m_changes = _changed_map[std::type_index(type(mgr))];
        auto &changed = m_changes[mgr];
        for (const auto &p: ptrs)
            changed.insert(p);
    }

    void clear()
    {
        _created_map.clear();
        _changed_map.clear();
    }

    

private:
    std::unordered_map<int, std::string> _reason_strings;
    ptr_type_to_set _created_map;
    ptr_type_to_set _changed_map;
    std::set<int> _reasons;
}; //class Change_Tracker
} //namespace isolde

#endif //ISOLDE_CHANGETRACKER
