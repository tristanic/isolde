#ifndef CLIPPER_SYMMETRY
#define CLIPPER_SYMMETRY

#include <iostream>

#include <cstddef>
#include <vector>
#include <array>
#include <unordered_map>
#include <utility>
#include <cmath>

template <typename T>
inline void transform_coord (T tf[3][4], T coord[3], T out[3])
{
    for (size_t i=0; i<3; ++i) {
        T* row = &tf[i];
        out[i] = row[0]*coord[0] + row[1]*coord[1] + row[2]*coord[2] + row[3];
    }
}

template <typename T>
inline void transform_coord (T tf[12], T coord[3], T out[3])
{
    for (size_t i=0; i<3; ++i) {
        T* row = tf + i*4;
        out[i] = row[0]*coord[0] + row[1]*coord[1] + row[2]*coord[2] + row[3];
    }
}

template <typename T>
inline T l2_norm_3d(T a[3])
{
    T accum = 0;
    for (int i = 0; i < 3; i++) {
        accum += a[i]*a[i];
    }
    return sqrt(accum);
}

template <typename T>
inline bool distance_below_cutoff(T a[3], T b[3], T cutoff)
{
    T deltas[3];
    for (size_t i=0; i<3; ++i) {
        deltas[i] = std::abs((*a++)-(*b++));
        if (deltas[i] > cutoff)
            return false;
    }
    if (l2_norm_3d(deltas) > cutoff)
        return false;
    return true;
}

template <typename T>
inline T distance_3D(T a[3], T b[3])
{
    T accum = 0;
    for (size_t i=0; i<3; ++i) {
        T diff = *(a++) - *(b++);
        accum += diff*diff;
    }
    return sqrt(accum);
}

template <typename T>
inline void copy_coord(T from[3], T to[3]) {
    for (size_t i=0; i<3; ++i)
        *(to++) = *(from++);
}

class Sym_Close_Points
{
public:
    typedef std::array<double, 3> Coord;
    typedef std::vector<size_t> indices;
    typedef std::vector<Coord> coords;
    typedef std::unordered_map<size_t, std::pair<indices, coords>> Sym_Map;
    Sym_Map symmap;

}; // class Sym_Closest_Points

const size_t TF_SIZE = 3*4;
void find_close_points_sym(double center[3], double cutoff, double *transforms, size_t n_sym,
    double *point_coords, size_t n_points, Sym_Close_Points &close)
{
    double tf_coord[3];
    for (size_t i=0; i<n_sym; ++i)
    {
        double *tf = transforms + i*TF_SIZE;
        double *pc = point_coords;
        Sym_Close_Points::indices close_indices;
        Sym_Close_Points::coords close_coords;
        for (size_t j=0; j<n_points; ++j) {
            transform_coord(tf, pc, tf_coord);

            if (distance_below_cutoff(tf_coord, center, cutoff))
            {
                Sym_Close_Points::Coord this_coord;
                copy_coord(tf_coord, this_coord.data());
                close_indices.push_back(j);
                close_coords.push_back(this_coord);
            }
            pc += 3;
        }
        close.symmap[i] = std::make_pair(close_indices, close_coords);
    }
}

#endif //CLIPPER_SYMMETRY
