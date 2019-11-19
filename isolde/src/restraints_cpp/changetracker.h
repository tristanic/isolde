/**
 * @Author: Tristan Croll <tic20>
 * @Date:   26-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 28-Mar-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#ifndef ISOLDE_CHANGETRACKER
#define ISOLDE_CHANGETRACKER

#include <vector>
#include <string>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <set>
#include <iostream>

#include <pyinstance/PythonInstance.declare.h>

namespace isolde
{
class PositionRestraintMgr;
class PositionRestraint;
class DistanceRestraintMgr;
class DistanceRestraint;
class DihedralRestraintMgr;
class Dihedral_Restraint;


class Change_Tracker: public pyinstance::PythonInstance<Change_Tracker>
{
public:
    /* Hierarchy:
     * <Manager Type>
     *      |
     *      ---<Manager ptr>
     *            |
     *            --- <Reason>
     *                   |
     *                   ----<Changed pointers>
     *
     */
    typedef std::unordered_map<int, std::set<const void*>> Reason_to_Ptrs;
    typedef std::unordered_map<const void*, Reason_to_Ptrs> Mgr_to_Managed;
    typedef std::unordered_map<std::type_index, Mgr_to_Managed> Ptr_Type_to_Set;

    enum Reasons{
        REASON_RESTRAINT_CREATED,
        REASON_TARGET_CHANGED,
        REASON_CUTOFF_CHANGED,
        REASON_SPRING_CONSTANT_CHANGED,
        REASON_ADAPTIVE_C_CHANGED,
        REASON_DISPLAY_CHANGED,
        REASON_ENABLED_CHANGED,
    };

    Change_Tracker();
    ~Change_Tracker() {}

    template<class Mgr>
    std::set<const void*>& get_change_set(std::type_index mgr_type, Mgr* mgr, int reason)
    {
        if (_python_class_name.find(mgr_type) == _python_class_name.end())
            throw std::out_of_range("Restraint manager type has not been registered!");

        return _change_map[mgr_type][static_cast<const void *>(mgr)][reason];
    }

    void register_mgr(const std::type_index &mgr_type,
        const std::string &mgr_pyname, const std::string &c_pyname)
    {
        _python_class_name[mgr_type] = std::make_pair(mgr_pyname, c_pyname);
    }

    const std::pair<std::string, std::string>& get_python_class_names(const std::type_index &mgr_type) const
    {
        return _python_class_name.at(mgr_type);
    }

    const std::string& reason_string(int reason) const
    {
        return _reason_strings.at(reason);
    }

    const std::unordered_map<int, std::string>& all_reason_strings() const
    {
        return _reason_strings;
    }

    template <class Mgr, class C>
    void add_created(const std::type_index &mgr_type, Mgr* mgr, C* ptr)
    {
        auto& cset = get_change_set(mgr_type, mgr, REASON_RESTRAINT_CREATED);
        cset.insert(ptr);
    }

    template <class Mgr, class C>
    void add_created_batch(const std::type_index &mgr_type, Mgr* mgr, std::vector<C*> ptrs)
    {
        auto& cset = get_change_set(mgr_type, mgr, REASON_RESTRAINT_CREATED);
        for (const auto &p: ptrs)
            cset.insert(p);
    }

    template <class Mgr, class C>
    void add_modified(const std::type_index &mgr_type, Mgr *mgr, C* ptr, int reason)
    {
        auto& cset = get_change_set(mgr_type, mgr, reason);
        cset.insert(ptr);
    }

    template <class Mgr, class C>
    void add_modified_batch(const std::type_index &mgr_type, Mgr *mgr, std::vector<C*> ptrs, int reason)
    {
        auto& cset = get_change_set(mgr_type, mgr, reason);
        for (const auto &p: ptrs)
            cset.insert(p);
    }

    void clear()
    {
        _change_map.clear();
    }


    const Ptr_Type_to_Set& get_changes() const {
        return _change_map;
    }




private:
    std::unordered_map<std::type_index, std::pair<std::string, std::string>> _python_class_name;
    std::unordered_map<int, std::string> _reason_strings;
    Ptr_Type_to_Set _change_map;

}; //class Change_Tracker
} //namespace isolde

#endif //ISOLDE_CHANGETRACKER
