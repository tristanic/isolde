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
    struct cutoffs 
    {
        double allowed;
        double log_allowed; 
        double outlier;
        double log_outlier;
        cutoffs() {}
        cutoffs(double a, double o): allowed(a), log_allowed(log(a)), outlier(o), log_outlier(log(o)) {}
    };
    //! RGBA colours in range 0..1
    typedef double color[4];
    struct colormap
    {
        color max; // most favoured
        color mid; // allowed cutoff
        color min; // least favoured
        color na; // Not scorable
    };
    
    
    void set_cutoffs(size_t r_case, const double &allowed, const double &outlier) {
        _cutoffs[r_case] = cutoffs(allowed, outlier);
    }
    cutoffs* get_cutoffs(size_t r_case) {return &(_cutoffs.at(r_case));}
    
    void set_colors(uint8_t *max, uint8_t *mid, uint8_t *min, uint8_t *na);
    colormap *get_colors();
    
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
    
    //! Score a single residue
    double validate(Residue *residue, Proper_Dihedral_Mgr *dmgr);
    
    void color_by_scores(double *scores, uint8_t *r_case, const size_t &n, uint8_t *colors);
    
private:
    std::unordered_map<size_t, RegularGridInterpolator<double>> _interpolators;
    std::unordered_map<size_t, cutoffs> _cutoffs;
    enum Rama_Case{NONE=0, CISPRO=1, TRANSPRO=2, GLYCINE=3, PREPRO=4, ILEVAL=5, GENERAL=6, NUM_RAMA_CASES=7};
    double _cis_cutoff = M_PI/6.0; 
    static const std::string OMEGA_STR;
    static const std::string PHI_STR;
    static const std::string PSI_STR;
    static const double NONE_VAL;
    colormap _colors;
    void _interpolate_colors(const color& min_color, const color& max_color, 
        const double &min_val, const double &max_val, const double &score, color &out);
    
        
}; //class Rama_Mgr    
}//namespace isolde

#endif //ISOLDE_RAMA
