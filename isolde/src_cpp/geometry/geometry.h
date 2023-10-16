/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 02-Apr-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#ifndef isolde_geometry
#define isolde_geometry

#include <math.h>
#include <cstdlib>
#include <cstdint>

#include "../constants.h"

namespace isolde
{
namespace geometry
{

template <typename T>
inline void cross_product_3D(T *a, T *b, T* result)
{
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}

template <typename T>
inline T dot_product_3D(T a[3], T b[3])
{
    T accum = 0;
    for (int i=0; i < 3; ++i) {
        accum += (*a++)*(*b++);
    }
    return accum;
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
inline void normalize_vector_3d(T vector[3])
{
    T norm = l2_norm_3d<T>(vector);
    for (int i = 0; i < 3; ++i) {
        vector[i] /= norm;
    }
}

// For direct use with coordinates
template <typename T>
inline T dihedral_angle(T p0[3], T p1[3], T p2[3], T p3[3])
{
    T b0[3];
    T b1[3];
    T b2[3];

    for (int i = 0; i < 3; i++) {
        b0[i] = p0[i] - p1[i];
        b1[i] = p2[i] - p1[i];
        b2[i] = p3[i] - p2[i];
    }

    T nb1 = l2_norm_3d(b1);

    for (int i =0; i < 3; i++) {
        b1[i] /= nb1;
    }

    T dp01 = dot_product_3D<T>(b0, b1);
    T dp21 = dot_product_3D<T>(b2, b1);

    T v[3];
    T w[3];

    for (int i = 0; i< 3; i++) {
        v[i] = b0[i] - dp01*b1[i];
        w[i] = b2[i] - dp21*b1[i];
    }

    T x = dot_product_3D<T>(v, w);
    T xp[3];
    cross_product_3D<T>(b1, v, xp);
    T y = dot_product_3D<T>(xp, w);

    return atan2(y, x);
}


// For use with ChimeraX Point objects
template <typename T, typename R>
inline R dihedral_angle(const T& p0, const T& p1, const T& p2, const T& p3)
{
    R b0[3];
    R b1[3];
    R b2[3];

    for (int i = 0; i < 3; i++) {
        b0[i] = p0[i] - p1[i];
        b1[i] = p2[i] - p1[i];
        b2[i] = p3[i] - p2[i];
    }

    R nb1 = l2_norm_3d(b1);

    for (int i =0; i < 3; i++) {
        b1[i] /= nb1;
    }

    R dp01 = dot_product_3D<R>(b0, b1);
    R dp21 = dot_product_3D<R>(b2, b1);

    R v[3];
    R w[3];

    for (int i = 0; i< 3; i++) {
        v[i] = b0[i] - dp01*b1[i];
        w[i] = b2[i] - dp21*b1[i];
    }

    R x = dot_product_3D<R>(v, w);
    R xp[3];
    cross_product_3D<R>(b1, v, xp);
    R y = dot_product_3D<R>(xp, w);

    return atan2(y, x);
}

template<typename T>
void multiply_transforms(T tf1[3][4], T tf2[3][4], T out[3][4])
{
    const T ll[4] = {0, 0, 0, 1};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            out[i][j] = tf1[i][0]*tf2[0][j] + tf1[i][1]*tf2[1][j]
                      + tf1[i][2]*tf2[2][j] + tf1[i][3]*ll[j];
        }
    }
} // multiply_transforms

template<typename T>
void multiply_transforms_gl(T* tf1, T* tf2, T* out)
{
    for (size_t i=0; i<4; ++i) {
        for(size_t j=0; j<4; ++j) {
            size_t row = 4*i;
            out[j+row] = tf1[row]*tf2[j] + tf1[row+1]*tf2[4+j]
                + tf1[row+2]*tf2[8+j] + tf1[row+3]*tf2[12+j];
        }
    }
}


// Fill out with a 3x4 matrix of rotation matrices (translation = 0)
template <typename T>
void rotation(T normalized_axis[3], const T &angle, T* out)
{
    T sa, ca, k, ax, ay, az;
    ax = normalized_axis[0]; ay = normalized_axis[1]; az = normalized_axis[2];
    sa = sin(angle);
    ca = cos(angle);
    k = 1 - ca;

    *out++ = 1 + k*( ax*ax -1);
    *out++ = -az*sa + k*ax*ay;
    *out++ = ay*sa + k*ax*ax;
    *out++ = 0;
    *out++ = az*sa + k*ax*ay;
    *out++ = 1 + k*(ay*ay - 1);
    *out++ = -ax*sa + k*ay*az;
    *out++ = 0;
    *out++ = -ay*sa + k*ax*az;
    *out++ = ax*sa + k*ay*az;
    *out++ = 1 + k*(az*az - 1);
    *out++ = 0;

} // rotations

//! Fill out with a 4x4 OpenGL rotation matrix about the given axis
/*! NOTE: expects the rotation axis to be already normalised!
*/
template <typename T>
void rotation_gl(const T normalized_axis[3], const T &angle, T* out)
{
    const T ca = cos(angle);
    const T sa = sin(angle);
    T ax = normalized_axis[0], ay = normalized_axis[1], az = normalized_axis[2];
    const T k = 1-ca;
    *out++ = 1 + k*( ax*ax -1);
    *out++ = az*sa + k*ax*ay;
    *out++ = -ay*sa + k*ax*az;
    *out++=0;

    *out++ = -az*sa + k*ax*ay;
    *out++ = 1 + k*(ay*ay - 1);
    *out++ = ax*sa + k*ay*az;
    *out++=0;

    *out++ = ay*sa + k*ax*ax;
    *out++ = -ax*sa + k*ay*az;
    *out++ = 1 + k*(az*az - 1);
    *out++ = 0;

    *out++ = 0;
    *out++ = 0;
    *out++ = 0;
    *out++ = 1;
}

//! Flip a 4x4 OpenGL transformation matrix on the x axis.
template <typename T>
void flip_on_x_gl(T *mat44)
{
    mat44[1]*=-1;
    mat44[2]*=-1;
    mat44[5]*=-1;
    mat44[6]*=-1;
    mat44[9]*=-1;
    mat44[10]*=-1;
    mat44[13]*=-1;
    mat44[14]*=-1;
}


template <class C, typename T>
void bond_cylinder_transform_gl(C xyz0, C xyz1, T r, T length_scale, T *rot44)
{
    T bvec[3];
    for (size_t i=0; i<3; i++) {
        bvec[i] = xyz1[i]-xyz0[i];
    }
    T d = l2_norm_3d(bvec);
    if (d < ALMOST_ZERO) {
        bvec[0]=0; bvec[1]=0; bvec[2]=1;
        d = ALMOST_ZERO;
    } else {
        for (auto &b: bvec)
            b/=d;
    }
    T &c = bvec[2], c1;
    if (c <= -1) c1 = 0;
    else c1 = 1.0/(1.0+c);
    T wx = -bvec[1], wy = bvec[0];
    T cx = c1*wx, cy = c1*wy;
    T h = d * length_scale;
    *rot44++ = r*(cx*wx + c);
    *rot44++ = r*cy*wx;
    *rot44++ = -r*wy;
    *rot44++ = 0;

    *rot44++ = r*cx*wy;
    *rot44++ = r*(cy*wy + c);
    *rot44++ = r*wx;
    *rot44++ = 0;

    *rot44++ = h*wy;
    *rot44++ = -h*wx;
    *rot44++ = h*c;
    *rot44++ = 0;

    *rot44++ = (xyz0[0]+xyz1[0])/2;
    *rot44++ = (xyz0[1]+xyz1[1])/2;
    *rot44++ = (xyz0[2]+xyz1[2])/2;
    *rot44++ = 1;
}



} // namespace geometry
} // namespace isolde
#endif //isolde_geometry
