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
class Position_Restraint_Mgr;
class Position_Restraint;
class Distance_Restraint_Mgr;
class Distance_Restraint;
class Dihedral_Restraint_Mgr;
class Dihedral_Restraint;


class Change_Tracker: public pyinstance::PythonInstance<Change_Tracker>
{
public:
    typedef std::unordered_map<const void*, std::set<const void *>> Mgr_To_Managed;
    typedef std::unordered_map<std::type_index, Mgr_To_Managed> Ptr_Type_To_Set;

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

    void register_mgr(const std::type_index &mgr_type,
        const std::string &mgr_pyname, const std::string &c_pyname)
    {
        _python_class_name[mgr_type] = std::make_pair(mgr_pyname, c_pyname);
    }

    template <class Mgr, class C>
    void add_created(Mgr* mgr, C* ptr)
    {
        auto &m_created = _created_map[std::type_index(typeid(mgr))];
        m_created[mgr].insert(ptr);
    }

    template <class Mgr, class C>
    void add_created_batch(Mgr* mgr, std::vector<C*> ptrs)
    {
        auto &m_created = _created_map[std::type_index(typeid(mgr))];
        auto &created = m_created[mgr];
        for (const auto &p: ptrs)
            created.insert(p);
    }

    template <class Mgr, class C>
    void add_modified(Mgr *mgr, C* ptr)
    {
        auto &m_changes = _changed_map[std::type_index(typeid(mgr))];
        m_changes[mgr].insert(ptr);
    }

    template <class Mgr, class C>
    void add_modified_batch(Mgr *mgr, std::vector<C*> ptrs)
    {
        auto &m_changes = _changed_map[std::type_index(typeid(mgr))];
        auto &changed = m_changes[mgr];
        for (const auto &p: ptrs)
            changed.insert(p);
    }

    void clear()
    {
        _created_map.clear();
        _changed_map.clear();
    }

    const Ptr_Type_To_Set& get_created() const {
        return _created_map;
    }

    const Ptr_Type_To_Set& get_changed() const {
        return _changed_map;
    }




private:
    std::unordered_map<std::type_index, std::pair<std::string, std::string>> _python_class_name;
    std::unordered_map<int, std::string> _reason_strings;
    Ptr_Type_To_Set _created_map;
    Ptr_Type_To_Set _changed_map;
    std::set<int> _reasons;
}; //class Change_Tracker
} //namespace isolde

#endif //ISOLDE_CHANGETRACKER
