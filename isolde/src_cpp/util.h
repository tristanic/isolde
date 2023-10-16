/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#ifndef ISOLDE_UTIL
#define ISOLDE_UTIL

#include <atomstruct/Coord.h>
#include "constants.h"

namespace isolde
{
namespace util
{
    inline void copy_coord(const atomstruct::Coord& coord_from, atomstruct::Coord coord_to)
    {
        coord_to[0] = coord_from[0];
        coord_to[1] = coord_from[1];
        coord_to[2] = coord_from[2];

    }


    inline double wrapped_angle(const double &angle) { return remainder(angle, TWO_PI); }
    inline double wrapped_angle_deg(const double &angle) { return remainder(angle, 360.0); }
    inline double degrees(const double &angle_rad) { return angle_rad * 180.0 / M_PI; }
    inline double radians(const double &angle_deg) { return angle_deg * M_PI / 180.0; }


} //namespace util
} //namespace isolde

#endif //ISOLDE_UTIL
