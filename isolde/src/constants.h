/**
 * @Author: Tristan Croll
 * @Date:   06-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
 * @Copyright: Copyright 2017-2018 Tristan Croll
 */



#ifndef ISOLDE_CONSTANTS
#define ISOLDE_CONSTANTS

#include <cmath>
namespace isolde
{

const double ALMOST_ZERO = 1e-12;
const double NONE_VAL = std::nan("Not applicable");
const double NO_RAMA_SCORE = -1.0;
const double TWO_PI = 2.0*M_PI;
const double NAN_NOT_SET = std::nan("Not set");
const int HIDE_ISOLDE = 0x02;
const double CIS_CUTOFF = M_PI/6.0;
const double TWISTED_CUTOFF = M_PI*5.0/6.0;
const double LINEAR_RESTRAINT_MAX_RADIUS = 0.3;
const double LINEAR_RESTRAINT_MIN_RADIUS = 0.025;
const double MAX_LINEAR_SPRING_CONSTANT = 100000.0;
const double MAX_RADIAL_SPRING_CONSTANT = 500.0;
const double DEFAULT_CHIRAL_RESTRAINT_SPRING_CONSTANT = 100.0;
const double DEFAULT_CHIRAL_RESTRAINT_CUTOFF = 15.0/180*M_PI;
const double DIHEDRAL_RESTRAINT_MAX_WIDTH = 3.0;
const double DIHEDRAL_RESTRAINT_MIN_WIDTH = 0.5;
const double Z_AXIS[3] = {0.0, 0.0, 1.0};
const double MIN_DISTANCE_RESTRAINT_TARGET = 1.0; //Angstroms

}

#endif //ISOLDE_CONSTANTS
