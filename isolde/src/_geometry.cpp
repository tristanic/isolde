#include <math.h>
#include <stdint.h>

#include "geometry/geometry.h"

typedef uint8_t npy_bool;
typedef float float32_t;
typedef double float64_t;


using namespace isolde::geometry;

extern "C"

{

double get_dihedral(double p0[3], double p1[3], double p2[3], double p3[3])
{
    return dihedral_angle<double>(p0, p1, p2, p3);
} // get_dihedral

void get_dihedrals(double *coords, int n, double * out)
{
    for (int i=0; i<n; ++i) {
        int j = i*12;
        *out++ = dihedral_angle<double>(coords+j, coords+j+3, coords+j+6, coords+j+9);
    }
} // get_dihedrals

// Fill out with a nx3x4 matrix of rotation matrices (translation = 0)
void rotations(double axis[3], double* angles, int n, double* out)
{
    normalize_vector_3d<double>(axis);
    for (int i = 0; i < n; ++i) {
        rotation<double>(axis, *angles++, out);
        out += 12;
    }
} // rotations

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
} // scale_transforms

void multiply_transforms(double tf1[3][4], double tf2[3][4], double out[3][4])
{
    multiply_transforms<double>(tf1, tf2, out);
}

// void multiply_transforms(double tf1[3][4], double tf2[3][4], double out[3][4])
// {
//     const double ll[4] = {0, 0, 0, 1};
//     for (int i = 0; i < 3; ++i) {
//         for (int j = 0; j < 4; ++j) {
//             out[i][j] = tf1[i][0]*tf2[0][j] + tf1[i][1]*tf2[1][j]
//                       + tf1[i][2]*tf2[2][j] + tf1[i][3]*ll[j];
//         }
//     }
// } // multiply_transforms

void flip_rotate_and_shift(int n, npy_bool* flip, double flip_op[3][4], double* rot, double* shift, double* out)
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
} // flip_rotate_and_shift

const int32_t T_INDICES[9] = {0,1,4,1,2,4,2,3,4};
const int D_LENGTH = 12;
const int V_LENGTH = 15;
const int N_LENGTH = 15;
const int T_LENGTH = 9;

// Generate vertices, normals and triangles to fill in a single dihedral
// with a pseudo-trapezoid.
void dihedral_fill_plane(double* coords, float* vertices, float* normals, int32_t* triangles, int start=0)
{
    int count = 0;
    for (int i=0; i<4; ++i) {
        for (int j=0; j<3; ++j, ++count) {
            vertices[count] = coords[count];
        }
    }
    for (int i = 0; i<3; ++i) {
        vertices[count+i] = (coords[i] + coords[9+i])/2;
    }

    // Pretend the surface is planar, and assign a single normal.
    float v1[3], v2[3];
    for (int i=0; i<3; ++i) {
        v1[i] = vertices[3+i]-vertices[i];
        v2[i] = vertices[9+i]-vertices[i];
    }

    cross_product_3D<float>(v1, v2, normals);

    // Copy the first normal over to the rest
    for (int i=3; i<N_LENGTH; ++i) {
        normals[i] = normals[i%3];
    }

    for (int i=0; i<9; ++i) {
        triangles[i] = start+T_INDICES[i];
    }
} // dihedral_fill_plane


void dihedral_fill_planes(int n, double* coords, float* vertices,
                            float* normals, int32_t* triangles)
{
    for (int i=0; i<n; ++i) {
        dihedral_fill_plane(coords+i*D_LENGTH,
                            vertices+i*V_LENGTH,
                            normals+i*N_LENGTH,
                            triangles+i*T_LENGTH,
                            i*V_LENGTH/3);
    }
} // dihedral_fill_planes

void dihedral_fill_and_color_planes(int n, double* coords, npy_bool* twisted_mask,
                                    npy_bool* cispro_mask, uint8_t* default_color,
                                    uint8_t* twisted_color, uint8_t* cispro_color,
                                    float* vertices, float* normals, int32_t* triangles,
                                    uint8_t* colors)
{
    dihedral_fill_planes(n, coords, vertices, normals, triangles);
    uint8_t* ptr;
    for (int i=0, count=0; i<n; ++i) {
        if (twisted_mask[i]) ptr = twisted_color;
        else if (cispro_mask[i]) ptr = cispro_color;
        else ptr = default_color;
        for (int j=0; j<5; ++j) {
            for (int k=0; k<4; ++k, count++) {
                colors[count] = ptr[k];
            }
        }
    }
} // dihedral_fill_and_color_planes

} // extern "C"
