/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */



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
private:
    int _sim_index = -1;
}; //SIM_RESTRAINT_BASE


} //namespace isolde


#endif //ISOLDE_SIM_RESTRAINT_BASE
