#include <math.h>

//using namespace std;

extern "C"

{

double* cross_product_3D(double a[3], double b[3])
{
    static double result[3];
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
    return result;
}

double dot_product_3D(double a[3], double b[3])
{
    double accum = 0;
    for (int i = 0; i < 3; i++) {
        accum += a[i]*b[i];
    }
    return accum;

}

double l2_norm_3d(double a[3])
{
    double accum = 0;
    for (int i = 0; i < 3; i++) {
        accum += a[i]*a[i];
    }
    return sqrt(accum);


}


double get_dihedral(double p0[3], double p1[3], double p2[3], double p3[3])
{
    double b0[3];
    double b1[3];
    double b2[3];

    for (int i = 0; i < 3; i++) {
        b0[i] = p0[i] - p1[i];
        b1[i] = p2[i] - p1[i];
        b2[i] = p3[i] - p2[i];
    }

    double nb1 = l2_norm_3d(b1);

    for (int i =0; i < 3; i++) {
        b1[i] /= nb1;
    }

    double dp01 = dot_product_3D(b0, b1);
    double dp21 = dot_product_3D(b2, b1);

    double v[3];
    double w[3];

    for (int i = 0; i< 3; i++) {
        v[i] = b0[i] - dp01*b1[i];
        w[i] = b2[i] - dp21*b1[i];
    }

    double x = dot_product_3D(v, w);
    double y = dot_product_3D(cross_product_3D(b1, v), w);

    return atan2(y, x);
}


}

