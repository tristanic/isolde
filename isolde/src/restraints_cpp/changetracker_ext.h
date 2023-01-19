/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#ifndef ISOLDE_CHANGES_EXT
#define ISOLDE_CHANGES_EXT

#include "changetracker.h"
#include "../molc.h"

using namespace atomstruct;
using namespace isolde;

/*******************************************************
 *
 * Change_Tracker functions
 *
 *******************************************************/
SET_PYTHON_INSTANCE(change_tracker, Change_Tracker)
GET_PYTHON_INSTANCES(change_tracker, Change_Tracker)

extern "C" EXPORT void*
change_tracker_new()
{
    try {
        Change_Tracker *t = new Change_Tracker();
        return t;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT void
change_tracker_delete(void *tracker)
{
    Change_Tracker *t = static_cast<Change_Tracker *>(tracker);
    try {
        delete t;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
change_tracker_clear(void *tracker)
{
    Change_Tracker *t = static_cast<Change_Tracker *>(tracker);
    try {
        t->clear();
    } catch (...) {
        molc_error();
    }
}

static PyObject* changes_as_py_dict(Change_Tracker *t,
    const Change_Tracker::Ptr_Type_to_Set &changes)
{
    PyObject* changes_data = PyDict_New();
    for (const auto &it1: changes)
    {
        const auto& py_classnames = t->get_python_class_names(it1.first);
        PyObject *mgr_type_key = unicode_from_string(py_classnames.first);
        PyObject* mgr_type_dict = PyDict_New();
        for (const auto &it2: it1.second)
        {
            PyObject *mgr_key = PyLong_FromVoidPtr(const_cast<void*>(it2.first));
            PyObject *change_dict = PyDict_New();
            for (const auto &it3: it2.second)
            {
                PyObject *change_key = unicode_from_string(t->reason_string(it3.first));
                void **ptrs;
                PyObject *ptr_array = python_voidp_array(it3.second.size(), &ptrs);
                for (auto ptr: it3.second)
                    (*ptrs++) = const_cast<void*>(ptr);
                PyDict_SetItem(change_dict, change_key, ptr_array);
                Py_DECREF(change_key);
                Py_DECREF(ptr_array);
            }
            PyDict_SetItem(mgr_type_dict, mgr_key, change_dict);
            Py_DECREF(mgr_key);
            Py_DECREF(change_dict);
        }
        PyDict_SetItem(changes_data, mgr_type_key, mgr_type_dict);
        Py_DECREF(mgr_type_key);
        Py_DECREF(mgr_type_dict);
    }
    return changes_data;
}

extern "C" EXPORT PyObject*
change_tracker_changes(void *tracker)
{
    Change_Tracker *t = static_cast<Change_Tracker *>(tracker);
    try {
        return changes_as_py_dict(t, t->get_changes());
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT PyObject*
change_tracker_reason_names(void *tracker)
{
    Change_Tracker *t = static_cast<Change_Tracker *>(tracker);
    try {
        const auto &reasons = t->all_reason_strings();
        size_t size = reasons.size();
        PyObject* ret = PyTuple_New(size);
        size_t i=0;
        for (const auto &r: reasons) {
            PyTuple_SET_ITEM(ret, i++, unicode_from_string(r.second));
        }
        return ret;
    } catch (...) {
        molc_error();
        return 0;
    }
}


#endif
