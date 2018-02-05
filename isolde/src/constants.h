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

}

#endif //ISOLDE_CONSTANTS
