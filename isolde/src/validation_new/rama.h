#ifndef ISOLDE_RAMA
#define ISOLDE_RAMA

#include <cmath>

#include "../atomic_cpp/dihedral.h"
#include "../atomic_cpp/dihedral_mgr.h"
#include "../interpolation/nd_interp.h"
#include <atomstruct/destruct.h>
#include <atomstruct/string_types.h>
#include <atomstruct/Residue.h>
#include <pyinstance/PythonInstance.declare.h>

namespace isolde
{
class Rama_Mgr: public pyinstance::PythonInstance<Rama_Mgr>
{
public:
    Rama_Mgr() {}
    ~Rama_Mgr() {}
    void add_interpolator(size_t r_case, const size_t &dim, 
        uint32_t *n, double *min, double *max, double *data);

    uint8_t rama_case(Dihedral *omega, Dihedral *psi);
    uint8_t rama_case(Residue *res, Proper_Dihedral_Mgr *mgr);

    void sort_into_rama_cases(Residue **residue, Dihedral **omega,
                  Dihedral **psi, const size_t &n, size_t *cases);
    
    //! Score a pre-processed set of residues
    /*!
     * This function is designed for highest-speed re-scoring of a set
     * of residues, where the phi, psi and omega values for each residue
     * have already been found, and Ramachandran cases have already been
     * determined. 
     */ 
    void validate(Residue **residue, Dihedral **omega, Dihedral **phi, 
                  Dihedral **psi, uint8_t *r_case, const size_t &n, double *scores);
    
    double validate(Residue *residue, Proper_Dihedral_Mgr *dmgr);
    
private:
    std::unordered_map<size_t, RegularGridInterpolator<double>> _interpolators;
    enum Rama_Case{NONE=0, CISPRO=1, TRANSPRO=2, GLYCINE=3, PREPRO=4, ILEVAL=5, GENERAL=6, NUM_RAMA_CASES=7};
    double _cis_cutoff = M_PI/6.0; 
    static const std::string OMEGA_STR;
    static const std::string PHI_STR;
    static const std::string PSI_STR;
    static const double NONE_VAL;
    
        
}; //class Rama_Mgr    
}//namespace isolde

#endif //ISOLDE_RAMA
