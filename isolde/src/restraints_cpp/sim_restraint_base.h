#ifndef ISOLDE_SIM_RESTRAINT_BASE
#define ISOLDE_SIM_RESTRAINT_BASE

namespace isolde
{

class Sim_Restraint_Base
{
public:
    void clear_sim_index() { _sim_index = -1; }
    void set_sim_index(int index) { _sim_index = index; }
    int get_sim_index() const { return _sim_index; }
    void set_sim_update_needed(bool flag) { _sim_update_needed = flag; }
    bool sim_update_needed() const { return _sim_update_needed; }
private:
    int _sim_index = -1;
    bool _sim_update_needed = true;
}; //SIM_RESTRAINT_BASE


} //namespace isolde


#endif //ISOLDE_SIM_RESTRAINT_BASE
