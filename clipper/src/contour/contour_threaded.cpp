#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <memory>
#include <future>
#include <iostream>

#include "contour.h"

using namespace Contour_Calculation;
namespace py=pybind11;

struct Contour_Geometry
{
public:
    float* vertex_xyz;
    float* normals;
    int* tv_indices;
    int vertex_count=0;
    int triangle_count=0;
    Contour_Geometry() {}
    Contour_Geometry(int vc, int tc)
        : vertex_count(vc), triangle_count(tc)
    {
        vertex_xyz = new float[vc*3];
        normals = new float[vc*3];
        tv_indices = new int[tc*3];
    }
};


class Contour_Thread_Mgr
{
public:
    Contour_Thread_Mgr() {}
    void start_compute(py::array_t<float> data, float threshold,
        bool cap_faces=true, bool return_normals=false);
    bool ready() const { return ready_; }
    bool return_normals() const { return return_normals_; }
    Contour_Geometry get_result();
private:
    Stride stride_[3];
    Index size_[3];
    std::unique_ptr<float> data_;
    std::future<Contour_Geometry> geom_;
    bool working_=false;
    bool ready_=false;
    bool return_normals_=false;
    void copy_3d_array(py::array_t<float> source);
    Contour_Geometry contour_surface_thread_(float* data, float threshold, bool cap_faces);
};

void Contour_Thread_Mgr::start_compute(py::array_t<float> data, float threshold,
    bool cap_faces, bool return_normals)
{
    if (working_)
        throw std::runtime_error("Contour thread is already running!");
    // Make a copy of the data to ensure C-contiguity and prevent it going out
    // of scope
    copy_3d_array(data);
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
    data_ = std::unique_ptr<float>(new float[source.shape(0)*source.shape(1)*source.shape(2)]);
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
    Contour_Geometry geom(cptr->vertex_count(), cptr->triangle_count());
    cptr->geometry(geom.vertex_xyz, reinterpret_cast<Index *>(geom.tv_indices));
    if (return_normals_)
        cptr->normals(geom.normals);
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
    T* d = reinterpret_cast<T *>(data);
    delete[] d;
}


PYBIND11_MODULE(contour_thread, m) {
    m.doc() = "Threaded contouring implementation";


    py::class_<Contour_Thread_Mgr>(m, "Contour_Thread_Mgr")
        .def(py::init<>())
        .def("start_compute", &Contour_Thread_Mgr::start_compute,
            py::arg("data"), py::arg("threshold"), py::arg("cap_faces")=true,
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
                 geom.vertex_xyz,
                 py::capsule(geom.vertex_xyz, &delete_when_done<float>));
            py::array_t<int> ta({tc, 3}, // shape
                {3*4, 4}, // C-style contiguous strides for n*3 int32
                geom.tv_indices,
                py::capsule(geom.tv_indices, &delete_when_done<int>));
            if (!self.return_normals())
                return py::make_tuple(va, ta);
            py::array_t<float> na({vc, 3}, // shape
                {3*4, 4}, // C-style contiguous strides for n*3 float
                geom.normals,
                py::capsule(geom.normals, &delete_when_done<float>));
            return py::make_tuple(va, ta, na);
        })
        ;


}
