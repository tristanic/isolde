/**
 * @Author: Tristan Croll <tic20>
 * @Date:   11-Jun-2019
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 11-Jun-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2016-2019 Tristan Croll
 */



#// General utilities not available through the official ChimeraX API

#include "../molc.h"

using namespace atomstruct;

extern "C" EXPORT PyObject*
residue_bonded_neighbors(void *residue)
{
    Residue *r = static_cast<Residue *>(residue);
    try {
        std::set<Residue *> rset;
        for (auto a: r->atoms())
        {
            auto neighbors = a->neighbors();
            for (auto n: neighbors)
            {
                if (n->residue() != r)
                    rset.insert(n->residue());
            }
        }
        void **rptr;
        PyObject* ra = python_voidp_array(rset.size(), &rptr);
        for (auto rr: rset)
            *(rptr++) = rr;
        return ra;
    } catch (...) {
        molc_error();
        return 0;
    }
}
