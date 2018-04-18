/**
 * @Author: Tristan Croll
 * @Date:   04-Feb-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   Tristan Croll
 * @Last modified time: 18-Apr-2018
 * @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
 * @Copyright: Copyright 2017-2018 Tristan Croll
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
