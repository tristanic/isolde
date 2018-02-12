#ifndef ISOLDE_CONSTANTS
#define ISOLDE_CONSTANTS

#include <cmath>
namespace isolde
{

const double NONE_VAL = std::nan("Not applicable");
const double NO_RAMA_SCORE = -1.0;
const double TWO_PI = 2.0*M_PI;
const double NAN_NOT_SET = std::nan("Not set");
const int HIDE_ISOLDE = 0x02;
const double CIS_CUTOFF = M_PI/6.0;
const double RESTRAINT_MAX_RADIUS = 0.4;
const double RESTRAINT_MIN_RADIUS = 0.025;
const double MAX_SPRING_CONSTANT = 5000.0;
const double Z_AXIS[3] = {0.0, 0.0, 1.0};

}

#endif //ISOLDE_CONSTANTS
