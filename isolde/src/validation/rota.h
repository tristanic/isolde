/**
 * @Author: Tristan Croll
 * @Date:   04-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   Tristan Croll
 * @Last modified time: 18-Apr-2018
 * @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
 * @Copyright: Copyright 2017-2018 Tristan Croll
 */



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

struct Rota_Target
{
    std::string name; // e.g. "mmtt"
    double frequency; // fractional rate of occurrence in high-resolution structures
    std::vector<double> angles; // target chi angles in radians
    std::vector<double> esds; // estimated standard deviations on each chi angle
    Rota_Target() {}
    Rota_Target(const std::string& tname, size_t nchi, double freq, double* ang, double* esd)
        : name(tname), frequency(freq)
    {
        for(size_t i=0; i<nchi; ++i)
        {
            angles.push_back(*ang++);
            esds.push_back(*esd++);
        }
    }
    bool operator<(const Rota_Target& other) const { return frequency < other.frequency; }
};

class Rota_Def
{
public:
    Rota_Def() {}
    Rota_Def(size_t n, size_t v, bool sym, const std::vector<std::vector<std::string>>& moving_atom_names)
        : _n_chi(n), _val_nchi(v), _symmetric(sym), _moving_atom_names(moving_atom_names) {}
    void add_target(const std::string& name, double freq, double* angles, double* esds);

    size_t n_chi() const { return _n_chi; }
    size_t val_nchi() const { return _val_nchi; }
    bool symmetric() const { return _symmetric; }
    const std::vector<Rota_Target>& targets() const { return _targets; }
    const std::vector<std::string>& moving_atom_names(size_t i) const
    {
        if (i >= _n_chi)
            throw std::out_of_range("This rotamer does not have that many chi dihedrals!");
        return _moving_atom_names[i];
    }
    // Sort targets in descending order of probability
    void sort_targets() { std::sort(_targets.rbegin(), _targets.rend()); }
    size_t num_targets() const { return _targets.size(); }
    Rota_Target* get_target(size_t i)
    {
        if (i >= num_targets())
            throw std::out_of_range("This rotamer does not have that many targets!");
        return &(_targets[i]);
    }

private:
    size_t _n_chi;
    size_t _val_nchi;
    bool _symmetric;
    std::vector<Rota_Target> _targets;
    std::vector<std::vector<std::string>> _moving_atom_names;
};



class Rotamer:  public pyinstance::PythonInstance<Rotamer>
{

public:
    Rotamer() {} // null constructor
    ~Rotamer() { auto du = DestructionUser(this); }
    Rotamer(Residue *res, Rota_Mgr *mgr);

    const std::vector<Proper_Dihedral *> &dihedrals() {return _chi_dihedrals; }
    size_t n_chi() const { return _def->n_chi(); }
    void angles(std::vector<double> &angles) const;
    void angles(double *angles) const;
    std::vector<double> angles() const;
    double score() const;
    Residue* residue() const { return _residue; }
    Structure* structure() const { return residue()->structure(); }
    Bond* ca_cb_bond() const { return _chi_dihedrals[0]->axial_bond(); }
    bool is_symmetric() const { return _def->symmetric(); }
    bool visible() const { return ca_cb_bond()->shown(); }

    Rota_Def* def() const { return _def; }
    size_t num_target_defs() const { return _def->num_targets(); }
    Rota_Target* get_target_def(size_t i) const { return _def->get_target(i); }

private:
    Residue* _residue;
    Rota_Mgr* _mgr;
    std::vector<Proper_Dihedral *> _chi_dihedrals;
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

    void add_rotamer_def(const std::string &resname, size_t n_chi, size_t val_nchi,
        bool symmetric, const std::vector<std::vector<std::string>>& moving_atom_names);
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
