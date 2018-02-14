#ifndef ISOLDE_ROTA
#define ISOLDE_ROTA

#include <thread> //testing

#include <string>
#include "../atomic_cpp/dihedral.h"
#include "../atomic_cpp/dihedral_mgr.h"
#include "../interpolation/nd_interp.h"
#include "../colors.h"
#include "../geometry/geometry.h"
#include <atomstruct/destruct.h>
#include <atomstruct/string_types.h>
#include <pyinstance/PythonInstance.declare.h>
#include <atomstruct/Residue.h>
#include <atomstruct/Bond.h>
using namespace atomstruct;

namespace isolde
{

class Rota_Mgr;

struct Rota_Def
{
    size_t n_chi;
    bool symmetric;
    Rota_Def() {}
    Rota_Def(size_t n, bool sym): n_chi(n), symmetric(sym) {}
};


class Rotamer:  public pyinstance::PythonInstance<Rotamer>
{

public:
    Rotamer() {} // null constructor
    ~Rotamer() { auto du = DestructionUser(this); }
    Rotamer(Residue *res, Rota_Mgr *mgr);

    const std::vector<Dihedral *> &dihedrals() {return _chi_dihedrals; }
    const size_t& n_chi() const { return _def->n_chi; }
    void angles(std::vector<double> &angles) const;
    void angles(double *angles) const;
    std::vector<double> angles() const;
    double score() const;
    Residue* residue() const {return _residue;}
    Bond* ca_cb_bond() const { return _chi_dihedrals[0]->axial_bond(); }
    bool is_symmetric() const { return _def->symmetric; }
    bool visible() const { return ca_cb_bond()->shown(); }

private:
    Residue* _residue;
    Rota_Mgr* _mgr;
    std::vector<Dihedral *> _chi_dihedrals;
    Rota_Def *_def;
    // size_t _n_chi;
    // bool _symmetric = false;

}; // class Rotamer

class Rota_Mgr: public DestructionObserver, public pyinstance::PythonInstance<Rota_Mgr>
{

public:
    enum Rota_Bins {FAVORED=0, ALLOWED=1, OUTLIER=2, BIN_NA=-1};
    Rota_Mgr() {} // null constructor
    Rota_Mgr(Proper_Dihedral_Mgr *dmgr): _dmgr(dmgr) {};
    ~Rota_Mgr();

    struct cutoffs
    {
        double allowed;
        double log_allowed;
        double outlier;
        double log_outlier;
        cutoffs() {}
        cutoffs(double a, double o): allowed(a), log_allowed(log(a)), outlier(o), log_outlier(log(o)) {}
    };

    void set_cutoffs(const double &allowed, const double &outlier) {_cutoffs = cutoffs(allowed, outlier);}
    cutoffs* get_cutoffs() {return &_cutoffs;}

    void set_colors(uint8_t *max, uint8_t *mid, uint8_t *min);
    colors::colormap *get_colors() {return &_colors;}

    void add_rotamer_def(const std::string &resname, size_t n_chi, bool symmetric);
    Rota_Def* get_rotamer_def(const std::string &resname);
    Rota_Def* get_rotamer_def(const ResName &resname);
    Rotamer* new_rotamer(Residue* residue);
    Rotamer* get_rotamer(Residue* residue);

    void add_interpolator(const std::string &resname, const size_t &dim,
        uint32_t *n, double *min, double *max, double *data);
    RegularGridInterpolator<double>* get_interpolator(const std::string &resname)
    {
        return &(_interpolators.at(resname));
    }
    RegularGridInterpolator<double>* get_interpolator(const ResName &resname)
    {
        return &(_interpolators.at(std::string(resname)));
    }
    Proper_Dihedral_Mgr* dihedral_mgr() { return _dmgr; }
    void validate(Rotamer** rotamers, size_t n, double* scores);
    void validate(Residue** residues, size_t n, double* scores);

    /**********TESTING***********/
    void validate_threaded(Rotamer **rotamers, size_t n, double* scores);
    bool thread_running() const { return _thread_running; }
    bool thread_done() const { return _thread_done; }
    void finalize_thread() { _validation_thread.join(); _thread_running=false; }
    /******END TESTING***********/


    int32_t bin_score(const double &score);
    void color_by_score(double *score, size_t n, uint8_t *out);
    virtual void destructors_done(const std::set<void*>& destroyed);

private:
    Proper_Dihedral_Mgr* _dmgr;
    std::unordered_map<Residue*, Rotamer*> _residue_to_rotamer;
    std::unordered_map<std::string, Rota_Def> _resname_to_rota_def;
    std::unordered_map<std::string, RegularGridInterpolator<double>> _interpolators;
    colors::colormap _colors;
    cutoffs _cutoffs;

    /*************TESTING********/
    void _validate_from_thread(Rotamer **rotamers, size_t n, double* scores);
    std::thread _validation_thread;
    bool _thread_done = false;
    bool _thread_running = false;

    /*********END TESTING********/


}; // class Rota_Mgr
} //namespace isolde

#endif //ISOLDE_ROTA
