/**
 * @Author: Tristan Croll <tic20>
 * @Date:   25-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 06-Jan-2020
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2016-2019 Tristan Croll
 */



#ifndef DIHEDRAL_MGR_EXT
#define DIHEDRAL_MGR_EXT

#include "dihedral_mgr.h"

#include "../molc.h"
using namespace atomstruct;
using namespace isolde;

/*************************************
 *
 * ProperDihedralMgr functions
 *
 *************************************/

SET_PYTHON_INSTANCE(proper_dihedral_mgr, ProperDihedralMgr)
GET_PYTHON_INSTANCES(proper_dihedral_mgr, ProperDihedralMgr)

extern "C" EXPORT void*
proper_dihedral_mgr_new()
{
   try {
       ProperDihedralMgr *mgr = new ProperDihedralMgr();
       return mgr;

   } catch (...) {
       molc_error();
       return nullptr;
   }
}

extern "C" EXPORT void
proper_dihedral_mgr_delete(void *mgr)
{
   ProperDihedralMgr *m = static_cast<ProperDihedralMgr *>(mgr);
   try {
       delete m;
   } catch (...) {
       molc_error();
   }
}

extern "C" EXPORT void
proper_dihedral_mgr_delete_dihedral(void *mgr, size_t n, void *dihedrals)
{
   ProperDihedralMgr *m = static_cast<ProperDihedralMgr *>(mgr);
   ProperDihedral **d = static_cast<ProperDihedral **>(dihedrals);
   try {
       std::set<ProperDihedral *> delete_list;
       for (size_t i=0; i<n; ++i) {
           delete_list.insert(d[i]);
       }
       m->delete_dihedrals(delete_list);
   } catch (...) {
       molc_error();
   }
}

extern "C" EXPORT void
proper_dihedral_mgr_add_dihedral_def(void *mgr, pyobject_t *rname,
   pyobject_t *dname, pyobject_t *anames, npy_bool *externals)
{
   ProperDihedralMgr *m = static_cast<ProperDihedralMgr *>(mgr);
   try {
       std::string resname(PyUnicode_AsUTF8(static_cast<PyObject *>(rname[0])));
       std::string dihe_name(PyUnicode_AsUTF8(static_cast<PyObject *>(dname[0])));
       std::vector<std::string> atom_names;
       std::vector<bool> externals_bool;
       for (size_t i=0; i<4; ++i) {
           atom_names.push_back(std::string(PyUnicode_AsUTF8(static_cast<PyObject *>(anames[i]))));
           externals_bool.push_back(externals[i]);
       }
       m->add_dihedral_def(resname, dihe_name, atom_names, externals_bool);
   } catch (...) {
       molc_error();
   }
}

extern "C" EXPORT void
proper_dihedral_mgr_reserve_map(void *mgr, size_t n)
{
   ProperDihedralMgr *m = static_cast<ProperDihedralMgr *>(mgr);
   try {
       m->reserve(n);
   } catch (...) {
       molc_error();
   }
}

extern "C" EXPORT void
proper_dihedral_mgr_new_dihedral(void *mgr, void*residues, size_t n, pyobject_t *name)
{
   ProperDihedralMgr *m = static_cast<ProperDihedralMgr *>(mgr);
   Residue **r = static_cast<Residue **>(residues);
   try {
       std::string sname = std::string(PyUnicode_AsUTF8(static_cast<PyObject *>(name[0])));
       for (size_t i=0; i<n; ++i) {
           m->new_dihedral(*r++, sname);
       }
   } catch (...) {
       molc_error();
   }
}


extern "C" EXPORT PyObject*
proper_dihedral_mgr_get_dihedrals(void *mgr, void *residues, pyobject_t *name, size_t n, npy_bool create)
{
   ProperDihedralMgr *m = static_cast<ProperDihedralMgr *>(mgr);
   Residue **r = static_cast<Residue **>(residues);
   try {
       std::vector<ProperDihedral *> dvec;
       std::string dname(PyUnicode_AsUTF8(static_cast<PyObject *>(name[0])));
       for (size_t i=0; i<n; ++i) {
           try {
               ProperDihedral *d = m->get_dihedral(r[i], dname, (bool)create);
               if (d != nullptr)
                   dvec.push_back(d);
           } catch (std::out_of_range&) {continue;}
       }
       void **dptr;
       PyObject *da = python_voidp_array(dvec.size(), &dptr);
       for (auto d: dvec) {
           *(dptr++) = d;
       }
       return da;
   } catch (...) {
       molc_error();
       return 0;
   }
}

extern "C" EXPORT void
proper_dihedral_mgr_residue_has_dihedral(void *mgr, void *residues, pyobject_t *name, size_t n, npy_bool* has_dihedral)
{
    ProperDihedralMgr *m = static_cast<ProperDihedralMgr *>(mgr);
    Residue **r = static_cast<Residue **>(residues);
    try {
        std::string dname(PyUnicode_AsUTF8(static_cast<PyObject *>(name[0])));
        for (size_t i=0; i<n; ++i) {
            try {
                has_dihedral[i] = (m->get_dihedral(r[i], dname, true) != nullptr);
            } catch (std::out_of_range&) {
                has_dihedral[i] = false;
            }
        }
    } catch (...) {
        molc_error();
    }
}


extern "C" EXPORT PyObject*
proper_dihedral_mgr_get_residue_dihedrals(void *mgr, void *residues, size_t n)
{
   ProperDihedralMgr *m = static_cast<ProperDihedralMgr *>(mgr);
   Residue **r = static_cast<Residue **>(residues);
   try {
       std::vector<ProperDihedral *> dvec;
       for (size_t i=0; i<n; ++i) {
           auto r_dvec = m->get_dihedrals(*r++);
           for (auto d: r_dvec)
               dvec.push_back(d);
       }
       void **dptr;
       PyObject *da = python_voidp_array(dvec.size(), &dptr);
       for (auto d: dvec)
           *(dptr++) = d;
       return da;
   } catch (...) {
       molc_error();
       return 0;
   }
}



extern "C" EXPORT int
proper_dihedral_mgr_num_mapped_dihedrals(void *mgr)
{
   ProperDihedralMgr *m = static_cast<ProperDihedralMgr *>(mgr);
   try {
       return m->num_mapped_dihedrals();
   } catch (...) {
       molc_error();
       return 0;
   }
}



#endif
