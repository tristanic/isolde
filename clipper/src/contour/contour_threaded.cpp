/**
 * @Author: Tristan Croll <tic20>
 * @Date:   05-Sep-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 15-May-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */



#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <memory>
#include <future>
#include <iostream>

// #include "bindings/numpy_helper.h"

#include "contour.h"

using namespace Contour_Calculation;
namespace py=pybind11;

struct Contour_Geometry
{
public:
    std::vector<float>* vertex_xyz;
    std::vector<float>* normals;
    std::vector<int>* tv_indices;
    // float* vertex_xyz;
    // float* normals;
    // int* tv_indices;
    int vertex_count=0;
    int triangle_count=0;
    Contour_Geometry() {}
    Contour_Geometry(int vc, int tc)
        : vertex_count(vc), triangle_count(tc)
    {
        vertex_xyz = new std::vector<float>(vc*3);
        normals = new std::vector<float>(vc*3);
        tv_indices = new std::vector<int>(tc*3);
        // vertex_xyz = new float[vc*3];
        // normals = new float[vc*3];
        // tv_indices = new int[tc*3];
    }
};

void copy_transform(py::array_t<float> source, float dest[3][4])
{
    for (size_t i=0; i<3; ++i)
        for (size_t j=0; j<4; ++j)
            dest[i][j] = *(source.data(i,j));
}


/* HELPER FUNCTIONS */

template <typename T>
void transform_coord (T tf[3][4], T coord[3], T out[3])
{
    for (size_t i=0; i<3; ++i) {
        T* row = tf[i];
        out[i] = row[0]*coord[0] + row[1]*coord[1] + row[2]*coord[2] + row[3];
    }
}

// Transform coordinates in-place
template <typename T>
void transform_coords (T tf[3][4], T* coords, int n)
{
    T temp[3];
    for (int i=0; i<n; ++i)
    {
        transform_coord(tf, coords, temp);
        for (int j=0; j<3; ++j)
            *coords++ = temp[j];
    }
}

template <typename T>
T l2_norm_3d(T a[3])
{
    T accum = 0;
    for (int i = 0; i < 3; i++) {
        accum += a[i]*a[i];
    }
    return sqrt(accum);
}

// Normalize 3D vectors in-place
template <typename T>
void normalize_3d_vectors(T* vecs, int n)
{
    for (int i=0; i<n; ++i)
    {
        T norm = l2_norm_3d(vecs);
        // Avoid divide-by-zero
        norm = norm > 0 ? norm : 1;
        for (int j=0; j<3; ++j)
            *vecs++ /= norm;
    }
}

// ----------------------------------------------------------------------------
// Swap vertex 1 and 2 of each triangle.
//
void reverse_triangle_vertex_order(int* ta, int n)
{
  int s0=3, s1 = 1;
  // int s0 = triangles.stride(0), s1 = triangles.stride(1);
  for (int t = 0 ; t < n ; ++t)
    {
      int i1 = s0*t+s1, i2 = i1 + s1;
      int v1 = ta[i1], v2 = ta[i2];
      ta[i1] = v2;
      ta[i2] = v1;
    }
}

class Contour_Thread_Mgr
{
public:
    Contour_Thread_Mgr() {}
    void start_compute(py::array_t<float> data, float threshold, float det,
        py::array_t<float> vertex_transform, py::array_t<float> normal_transform,
        bool cap_faces=true, bool return_normals=false);
    bool ready() const { return ready_; }
    bool return_normals() const { return return_normals_; }
    Contour_Geometry get_result();
private:
    std::unique_ptr<float[]> data_;
    float threshold_;
    bool flip_triangles_; // if det < 0
    float vertex_transform_[3][4]; // Transform mapping vertices into model coordinates
    float normal_transform_[3][4]; // Transform mapping normals into model coordinates

    Stride stride_[3];
    Index size_[3];
    std::future<Contour_Geometry> geom_;
    bool working_=false;
    bool ready_=false;
    bool return_normals_=false;
    void copy_3d_array(py::array_t<float> source);
    Contour_Geometry contour_surface_thread_(float data[], float threshold, bool cap_faces);
};

void Contour_Thread_Mgr::start_compute(py::array_t<float> data, float threshold,
    float det, py::array_t<float> vertex_transform, py::array_t<float> normal_transform,
    bool cap_faces, bool return_normals)
{
    if (working_)
        throw std::runtime_error("Contour thread is already running!");
    // Make a copy of the data to ensure C-contiguity and prevent it going out
    // of scope
    threshold_=threshold;
    copy_3d_array(data);
    flip_triangles_ = det < 0;
    copy_transform(vertex_transform, vertex_transform_);
    copy_transform(normal_transform, normal_transform_);
    working_=true;
    return_normals_ = return_normals;
    ready_=false;
    try {
        geom_ = std::async(std::launch::async, &Contour_Thread_Mgr::contour_surface_thread_, this, data_.get(), threshold, cap_faces);
    } catch (...) {
        working_ = false;
        throw;
    }
}

void Contour_Thread_Mgr::copy_3d_array(py::array_t<float> source)
{
    data_ = std::unique_ptr<float[]>(new float[source.shape(0)*source.shape(1)*source.shape(2)]);
    int stride_size = sizeof(float);
    // stride and size are stored reversed since ChimeraX keeps volume data in
    // z,y,x order while contour code requires x,y,z
    // NOTE: Numpy array strides are per byte regardless of data type, whereas
    // surface() wants strides per float
    for (size_t i=0; i<3; ++i)
    {
        stride_[i] = source.strides(2-i)/stride_size;
        size_[i] = source.shape(2-i);
    }
    auto dptr = data_.get();
    auto s = source.unchecked<3>();
    for (ssize_t i=0; i < source.shape(0); i++)
        for (ssize_t j=0; j < source.shape(1); j++)
            for (ssize_t k=0; k < source.shape(2); k++)
                dptr[i+stride_[1]*j+stride_[0]*k] = s(i,j,k);

}

Contour_Geometry Contour_Thread_Mgr::contour_surface_thread_(float* data, float threshold, bool cap_faces)
{
    auto cptr = std::unique_ptr<Contour_Surface>(surface(data, size_, stride_, threshold, cap_faces));
    int vc = cptr->vertex_count();
    int tc = cptr->triangle_count();
    Contour_Geometry geom(vc, tc);
    cptr->geometry(geom.vertex_xyz->data(), reinterpret_cast<Index *>(geom.tv_indices->data()));
    if (return_normals_)
        cptr->normals(geom.normals->data());
    if (flip_triangles_)
        reverse_triangle_vertex_order(geom.tv_indices->data(), tc);

    // Transform coords and normals to model coordinates
    transform_coords(vertex_transform_, geom.vertex_xyz->data(), vc);
    if (return_normals_)
    {
        transform_coords(normal_transform_, geom.normals->data(), vc);
        // Set all normals to unit length;
        normalize_3d_vectors(geom.normals->data(), vc);

    }
    ready_=true;
    return geom;
}

Contour_Geometry
Contour_Thread_Mgr::get_result()
{
    if (!working_)
        throw std::runtime_error("No contour calculation has been started!");
    working_ = false;
    return geom_.get();
}

template<typename T>
void delete_when_done(void *data)
{
    std::vector<T>* d = reinterpret_cast<std::vector<T> *>(data);
    delete d;
}



PYBIND11_MODULE(contour_thread, m) {
    m.doc() = "Threaded contouring implementation";


    py::class_<Contour_Thread_Mgr>(m, "Contour_Thread_Mgr")
        .def(py::init<>())
        .def("start_compute", &Contour_Thread_Mgr::start_compute,
            py::arg("data"), py::arg("threshold"), py::arg("determinant"),
            py::arg("vertex_transform"), py::arg("normal_transform"),
            py::arg("cap_faces")=true,
            py::arg("return_normals")=false)
        .def("ready", &Contour_Thread_Mgr::ready)
        .def("get_result", [](Contour_Thread_Mgr& self)
        {
            // Takes ownership of vertex_xyz, tv_indices and normals. They will
            // be automatically deleted when they go out of scope in Python.
            auto geom = self.get_result();
            const auto& vc = geom.vertex_count;
            const auto& tc = geom.triangle_count;
            py::array_t<float> va({vc, 3}, // shape
                 {3*4, 4}, // C-style contiguous strides for n*3 float
                 geom.vertex_xyz->data(),
                 py::capsule(geom.vertex_xyz, &delete_when_done<float>));
            py::array_t<int> ta({tc, 3}, // shape
                {3*4, 4}, // C-style contiguous strides for n*3 int32
                geom.tv_indices->data(),
                py::capsule(geom.tv_indices, &delete_when_done<int>));
            if (!self.return_normals())
                return py::make_tuple(va, ta);
            py::array_t<float> na({vc, 3}, // shape
                {3*4, 4}, // C-style contiguous strides for n*3 float
                geom.normals->data(),
                py::capsule(geom.normals, &delete_when_done<float>));
            return py::make_tuple(va, ta, na);
        })
        ;


}
