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

void get_dihedrals(double * coords, int n, double * out)
{
    double p0[3];
    double p1[3];
    double p2[3];
    double p3[3];

    for (int i = 0; i < n; i++) {
        int ii = 12*i;
        for (int j = 0; j < 3; j++) {
            p0[j] = coords[ii+j];
            p1[j] = coords[ii+3+j];
            p2[j] = coords[ii+6+j];
            p3[j] = coords[ii+9+j];
        }
        out[i] = get_dihedral(p0, p1, p2, p3);
    }

}

void normalize_vector_3d(double vector[3])
{
    double norm = l2_norm_3d(vector);
    for (int i = 0; i < 3; ++i) {
        vector[i] /= norm;
    }
}

// Fill out with a nx3x4 matrix of rotation matrices (translation = 0)
void rotations(double axis[3], double* angles, int n, double* out)
{
    normalize_vector_3d(axis);
    double angle, sa, ca, k, ax, ay, az;
    int start;
    ax = axis[0]; ay = axis[1]; az = axis[2];
    for (int i = 0; i < n; ++i) {
        start = i*12;
        angle = angles[i];
        sa = sin(angle);
        ca = cos(angle);
        k = 1 - ca;
        out[start  ] = 1 + k*( ax*ax -1);
        out[start+1] = -az*sa + k*ax*ay;
        out[start+2] = ay*sa + k*ax*ax;
        out[start+3] = 0;
        out[start+4] = az*sa + k*ax*ay;
        out[start+5] = 1 + k*(ay*ay - 1);
        out[start+6] = -ax*sa + k*ay*az;
        out[start+7] = 0;
        out[start+8] = -ay*sa + k*ax*az;
        out[start+9] = ax*sa + k*ay*az;
        out[start+10] = 1 + k*(az*az - 1);
        out[start+11] = 0;
    }
}

// Scale the input transforms by the values in scales, leaving the 
// translation components untouched.
void scale_transforms(double* scales, int n, double* transforms)
{
    for (int i = 0, count = 0; i < n; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 4; ++k, ++count)
            {
                if (k<3) transforms[count] *= scales[i];
            }
        }
    }
}


void multiply_transforms(double tf1[3][4], double tf2[3][4], double out[3][4])
{
    double ll[4] = {0, 0, 0, 1};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            out[i][j] = tf1[i][0]*tf2[0][j] + tf1[i][1]*tf2[1][j] 
                      + tf1[i][2]*tf2[2][j] + tf1[i][3]*ll[j];
        }
    }
    
}

void flip_rotate_and_shift(int n, bool* flip, double flip_op[3][4], double* rot, double* shift, double* out)
{
    double intermediate_result[3][4];
    double current_rot[3][4];
    double current_shift[3][4];
    double current_out[3][4];
    int tf_start;
    for (int i = 0; i < n; ++i) {
        tf_start = i*12;
        for (int j = 0, count = 0; j < 3; ++j) {
            for (int k = 0; k < 4; ++k, ++count) {
                current_rot[j][k] = rot[tf_start + count];
                current_shift[j][k] = shift[tf_start + count];
            }
        }
        
        if (flip[i]) {
            multiply_transforms(current_rot, flip_op, intermediate_result);
            multiply_transforms(current_shift, intermediate_result, current_out);
        } else {
            multiply_transforms(current_shift, current_rot, current_out);
        }
        for (int j = 0, count = 0; j < 3; ++j) {
            for (int k = 0; k < 4; ++k, ++count) {
                out[tf_start+count] = current_out[j][k];
            }
        }
        
        
    }
    
}

}

