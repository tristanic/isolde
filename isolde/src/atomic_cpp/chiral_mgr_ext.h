/**
 * @Author: Tristan Croll <tic20>
 * @Date:   25-Apr-2018
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 08-Nov-2020
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2016-2019 Tristan Croll
 */



#ifndef CHIRAL_MGR_EXT
#define CHIRAL_MGR_EXT

#include "chiral_mgr.h"

#include "../molc.h"
using namespace atomstruct;
using namespace isolde;

/*************************************
 *
 * ChiralMgr functions
 *
 *************************************/

SET_PYTHON_INSTANCE(chiral_mgr, ChiralMgr)
GET_PYTHON_INSTANCES(chiral_mgr, ChiralMgr)

extern "C" EXPORT void*
chiral_mgr_new()
{
   try {
       ChiralMgr *mgr = new ChiralMgr();
       return mgr;

   } catch (...) {
       molc_error();
       return nullptr;
   }
}

extern "C" EXPORT void
chiral_mgr_delete(void *mgr)
{
   ChiralMgr *m = static_cast<ChiralMgr *>(mgr);
   try {
       delete m;
   } catch (...) {
       molc_error();
   }
}

extern "C" EXPORT void
chiral_mgr_delete_chiral(void *mgr, size_t n, void *chirals)
{
   ChiralMgr *m = static_cast<ChiralMgr *>(mgr);
   ChiralCenter **c = static_cast<ChiralCenter **>(chirals);
   try {
       std::set<ChiralCenter *> delete_list;
       for (size_t i=0; i<n; ++i) {
           delete_list.insert(*c++);
       }
       m->delete_chirals(delete_list);
   } catch (...) {
       molc_error();
   }
}

/*
Chirals need to be treated with care, since unlike the standard named proper
dihedrals it is possible for a chiral centre to have different possible
substituent atoms when it forms an inter-residue vond. A good example is in
the glycosidic bonds of polysaccharides: almost every carbon is chiral, and the
name (and indeed element) of the bonded atom varies based on the bonding
pattern. The safest solution is to make each substituent in the definition be
a list of possible names.
*/

extern "C" EXPORT void
chiral_mgr_add_chiral_def(void *mgr, pyobject_t *rname, pyobject_t *cname,
    pyobject_t *s1_names, pyobject_t *s2_names, pyobject_t *s3_names,
    size_t ns1, size_t ns2, size_t ns3, double expected_angle, npy_bool* externals)
{
    ChiralMgr *m = static_cast<ChiralMgr *>(mgr);
    try {
        std::string resname(PyUnicode_AsUTF8(static_cast<PyObject *>(rname[0])));
        std::string chiral_name(PyUnicode_AsUTF8(static_cast<PyObject *>(cname[0])));
        std::vector<std::string> subs1;
        for (size_t i=0; i<ns1; ++i)
            subs1.push_back(PyUnicode_AsUTF8(static_cast<PyObject *>(*s1_names++)));
        std::vector<std::string> subs2;
        for (size_t i=0; i<ns2; ++i)
            subs2.push_back(PyUnicode_AsUTF8(static_cast<PyObject *>(*s2_names++)));
        std::vector<std::string> subs3;
        for (size_t i=0; i<ns3; ++i)
            subs3.push_back(PyUnicode_AsUTF8(static_cast<PyObject *>(*s3_names++)));
        std::vector<bool> ext;
        for (size_t i=0; i<3; ++i)
            ext.push_back(externals[i]);
        m->add_chiral_def(resname, chiral_name, subs1, subs2, subs3, expected_angle, ext);
    } catch(...) {
        molc_error();
    }
}

extern "C" EXPORT PyObject*
chiral_mgr_get_chiral(void *mgr, void *atom, size_t n, npy_bool create)
{
    ChiralMgr *m = static_cast<ChiralMgr *>(mgr);
    Atom **a = static_cast<Atom **>(atom);
    try {
        std::vector<ChiralCenter *> cvec;
        for (size_t i=0; i<n; ++i) {
            try {
                ChiralCenter *c = m->get_chiral(*a++, create);
                if (c!=nullptr)
                    cvec.push_back(c);
            } catch(std::out_of_range&) { continue; }
        }
        void **cptr;
        PyObject *ca = python_voidp_array(cvec.size(), &cptr);
        for (auto c: cvec)
            *(cptr++) = c;
        return ca;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT size_t
chiral_mgr_num_chirals(void *mgr)
{
    ChiralMgr *m = static_cast<ChiralMgr *>(mgr);
    try {
        return m->num_chirals();
    } catch (...) {
        molc_error();
        return 0;
    }
}

#endif //CHIRAL_MGR_EXT
