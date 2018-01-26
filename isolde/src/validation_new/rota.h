#ifndef ISOLDE_ROTA
#define ISOLDE_ROTA

#include "../atomic_cpp/dihedral.h"
#include "../atomic_cpp/dihedral_mgr.h"
#include "../interpolation/nd_interp.h"
#include <atomstruct/destruct.h>
#include <atomstruct/string_types.h>
#include <pyinstance/PythonInstance.declare.h>
#include <atomstruct/Residue.h>
#include <atomstruct/Bond.h>



using namespace atomstruct;

namespace isolde
{

class Rota_Mgr;

class Rotamer:  public pyinstance::PythonInstance<Rotamer>
{

public:
    Rotamer() {} // null constructor
    ~Rotamer() { auto du = DestructionUser(this); }
    Rotamer(Residue *res, Rota_Mgr *mgr);
    
    const std::vector<Dihedral *> &dihedrals() {return _chi_dihedrals; }
    const size_t& n_chi() const { return _n_chi; }
    void angles(std::vector<double> &angles);
    void angles(double *angles);
    std::vector<double> angles();
    double score();
    Residue* residue() {return _residue;}
    Bond* chi_axis() { return _chi_dihedrals[0]->axial_bond(); }
    
private:
    Residue* _residue;
    Rota_Mgr* _mgr;
    std::vector<Dihedral *> _chi_dihedrals;
    size_t _n_chi;

    
}; // class Rotamer

class Rota_Mgr: public DestructionObserver, public pyinstance::PythonInstance<Rota_Mgr>
{

public:
    Rota_Mgr() {} // null constructor
    ~Rota_Mgr();
    Rota_Mgr(Proper_Dihedral_Mgr *dmgr): _dmgr(dmgr) {};
    
    void add_rotamer_def(const std::string &resname, const std::vector<std::string> &chi_names);
    const std::vector<std::string>& get_rotamer_def(const std::string &resname);
    const std::vector<std::string>& get_rotamer_def(const ResName &resname);
    Rotamer* new_rotamer(Residue* residue);
    Rotamer* get_rotamer(Residue* residue);
    
    void add_interpolator(const std::string &resname, const size_t &dim, 
        uint32_t *n, double *min, double *max, double *data);
    const RegularGridInterpolator<double>& get_interpolator(const std::string &resname)
        { return _interpolators.at(resname); }
    const RegularGridInterpolator<double>& get_interpolator(const ResName &resname)
        {return _interpolators.at(std::string(resname)); }
    Proper_Dihedral_Mgr* dihedral_mgr() { return _dmgr; }
    void validate(Rotamer** rotamers, size_t n, double* scores);
    void validate(Residue** residues, size_t n, double* scores);
    
    virtual void destructors_done(const std::set<void*>& destroyed);

private:
    Proper_Dihedral_Mgr* _dmgr;
    std::unordered_map<Residue*, Rotamer*> _residue_to_rotamer;
    std::unordered_map<std::string, std::vector<std::string>> _resname_to_rota_def;
    std::unordered_map<std::string, RegularGridInterpolator<double>> _interpolators;

    
}; // class Rota_Mgr


} //namespace isolde

#endif //ISOLDE_ROTA
