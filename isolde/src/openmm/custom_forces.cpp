
#ifdef _WIN32
# define EXPORT __declspec(dllexport)
#else
# define EXPORT __attribute__((__visibility__("default")))
#endif

#include "../molc.h"
#include <vector>
#include <OpenMM.h>

extern "C"
{

EXPORT void
topoutbondforce_update_bond_parameters(void *force, size_t n, int *indices, double *k, double *r0)
{
    OpenMM::CustomBondForce *f = static_cast<OpenMM::CustomBondForce *>(force);
    try {
        std::vector<double> params(2);
        int particle1, particle2;
        for (size_t i =0; i<n; ++i) {
            int index = *(indices++);
            f->getBondParameters(index, particle1, particle2, params);
            params[0] = *k++;
            params[1] = *r0++;
            f->setBondParameters(index, particle1, particle2, params);
        }
    } catch (...) {
        molc_error();
    }
}

}
