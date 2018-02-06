
#ifndef isolde_geometry
#define isolde_geometry

#include <math.h>
#include <cstdint>

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
        accum += a[i]*b[i];
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


} // namespace geometry
} // namespace isolde
#endif //isolde_geometry
