/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#ifndef ISOLDE_RAMA
#define ISOLDE_RAMA

#include <cmath>
#include <array>

#include "../atomic_cpp/dihedral.h"
#include "../atomic_cpp/dihedral_mgr.h"
#include "../interpolation/nd_interp.h"
#include "../colors.h"

#include <atomstruct/destruct.h>
#include <atomstruct/string_types.h>
#include <atomstruct/Atom.h>
#include <atomstruct/Residue.h>
#include <pyinstance/PythonInstance.declare.h>

namespace isolde
{
class RamaMgr;

//! Ramachandran dihedrals for a single amino acid residue
class Rama: public pyinstance::PythonInstance<Rama>
{
public:
    Rama() {}
    ~Rama() { auto du = DestructionUser(this); }
    Rama(Residue* residue, ProperDihedralMgr *dmgr, RamaMgr *rmgr);
    //! Returns the omega dihedral, or nullptr if it doesn't exist
    Dihedral* omega();
    //! Returns the phi dihedral, or nullptr if it doesn't exist
    Dihedral* phi();
    //! Returns the psi dihedral, or nullptr if it doesn't exist
    Dihedral* psi();
    //! Returns the alpha carbon for the residue, or nullptr if it doesn't exist
    Atom* CA_atom() const { return _CA_atom; }
    //! Provides the omega, phi and psi angles, with std::nan for any that don't exist
    void angles(double *angles);
    void phipsi(double *angles);
    double omega_angle();
    double phi_angle();
    double psi_angle();
    double score();
    bool is_valid_rama();
    bool visible() const { return CA_atom()->visible(); }
    bool only_hidden_by_ribbon() const;
    uint8_t rama_case();
    bool no_valid_dihedrals();
    Residue *residue() const {return _residue;}
    bool check_for_deleted_dihedrals( const std::set<void *> &destroyed);

private:
    Dihedral* _omega = nullptr;
    Dihedral* _phi = nullptr;
    Dihedral* _psi = nullptr;
    Atom* _CA_atom;
    Residue* _residue;
    ProperDihedralMgr* _dmgr;
    RamaMgr* _rmgr;

    const char* err_msg_not_protein()
    {
        return "Residue must be part of a protein chain!";
    }
    const char* err_msg_no_dihedrals()
    {
        return "Residue must have the atoms to make at least one of phi, psi or omega!";
    }
    static const std::string OMEGA_STR;
    static const std::string PHI_STR;
    static const std::string PSI_STR;
    static const AtomName CA_NAME;

}; //Rama

//! Master manager for Ramachandran statistics and validation.
/*! Provides methods for getting/creating Rama objects for specific residues,
 *  scoring against MolProbity statistics, and colouring according to score.
 */
class RamaMgr: public pyinstance::PythonInstance<RamaMgr>, public DestructionObserver
{
public:
    RamaMgr() {}
    ~RamaMgr();
    RamaMgr(ProperDihedralMgr *mgr): _mgr(mgr) {}
    enum RamaCase{CASE_NONE=0, CISPRO=1, TRANSPRO=2, GLYCINE=3, PREPRO=4, ILEVAL=5, GENERAL=6, NUM_RAMA_CASES=7};
    enum RamaBins{FAVORED=0, ALLOWED=1, OUTLIER=2, BIN_NA=-1};
    struct cutoffs
    {
        double allowed;
        double log_allowed;
        double outlier;
        double log_outlier;
        cutoffs() {}
        cutoffs(double a, double o): allowed(a), log_allowed(log(a)), outlier(o), log_outlier(log(o)) {}
    };

    Rama* get_rama(Residue *res);

    void set_cutoffs(size_t r_case, const double &outlier, const double &allowed) {
        _cutoffs[r_case] = cutoffs(allowed, outlier);
    }
    cutoffs* get_cutoffs(size_t r_case) {return &(_cutoffs.at(r_case));}

    void set_colors(uint8_t *max, uint8_t *mid, uint8_t *min, uint8_t *na);
    colors::colormap *get_colors(size_t r_case) { return &(_colors.at(r_case)); }
    const colors::color& default_color() const {return _null_color;}

    void add_interpolator(size_t r_case, const size_t &dim,
        uint32_t *n, double *min, double *max, double *data);
    RegularGridInterpolator<double> *get_interpolator(size_t r_case)
        { return &(_interpolators.at(r_case)); }

    const colors::intcolor& cis_pro_color() const { return _cis_pro_color; }
    const colors::intcolor& cis_nonpro_color() const { return _cis_nonpro_color; }
    const colors::intcolor& twisted_color() const { return _twisted_color; }

    uint8_t rama_case(Residue *res);

    //! Get the Ramachandran P-value for a single residue
    double validate(Rama *r);
    //! Get the Ramachandran P-value for a single residue
    double validate(Residue *r);

    //! Get the Ramachandran scores for a set of residues
    /*! This could in principle be done by simply looping n times over
     *  validate(Rama*), but due to the way the RegularGridInterpolator
     *  is implemented it is more efficient to sort the residues into
     *  their cases first, then validate each set as a single array.
     */
    void validate(Rama **ramas, size_t n, double *scores, uint8_t *r_cases);

    //! Score a pre-processed set of residues
    /*!
     * This function is designed for highest-speed re-scoring of a set
     * of residues, where the phi, psi and omega values for each residue
     * have already been found, and Ramachandran cases have already been
     * determined.
     */

    void color_by_scores(double *scores, uint8_t *r_case, size_t n, uint8_t *out);
    int32_t bin_score(const double &score, uint8_t r_case);

    void delete_ramas(const std::set<Rama *> to_delete);
    virtual void destructors_done(const std::set<void*>& destroyed);

private:
    ProperDihedralMgr* _mgr;
    std::unordered_map<Residue*, Rama*> _residue_to_rama;
    std::unordered_map<size_t, RegularGridInterpolator<double>> _interpolators;
    std::unordered_map<size_t, cutoffs> _cutoffs;
    std::unordered_map<size_t, colors::colormap> _colors;
    colors::color _null_color;
    colors::intcolor _cis_pro_color = {64, 255, 64, 255};
    colors::intcolor _cis_nonpro_color = {255, 64, 64, 255};
    colors::intcolor _twisted_color = {255, 255, 64, 255};

    void _color_by_score(const double &score, const uint8_t &r_case, colors::color &color);
    void _delete_ramas(const std::set<Rama *> to_delete);
}; //class RamaMgr
}//namespace isolde

#endif //ISOLDE_RAMA
