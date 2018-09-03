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

py::array_t<float> copy_3d_array(py::array_t<float> source)
{
    py::array_t<float> ret({source.shape(0), source.shape(1), source.shape(2)});
    auto d = ret.mutable_unchecked<3>();
    auto s = source.unchecked<3>();
    for (ssize_t i=0; i < source.shape(0); i++)
        for (ssize_t j=0; j < source.shape(1); j++)
            for (ssize_t k=0; k < source.shape(2); k++)
                d(i,j,k) = s(i,j,k);
    return ret;
}

class Contour_Thread_Mgr
{
public:
    Contour_Thread_Mgr() {}
    void start_compute(py::array_t<float,0> data, float threshold,
        bool cap_faces=true, bool return_normals=false);
    bool ready() const { return ready_; }
    bool return_normals() const { return return_normals_; }
    Contour_Geometry get_result();
private:
    py::array_t<float> data_;
    std::future<Contour_Geometry> geom_;
    bool working_=false;
    bool ready_=false;
    bool return_normals_=false;
    Contour_Geometry contour_surface_thread_(py::array_t<float> data, float threshold, bool cap_faces);
};

void Contour_Thread_Mgr::start_compute(py::array_t<float,0> data, float threshold,
    bool cap_faces, bool return_normals)
{
    if (working_)
        throw std::runtime_error("Contour thread is already running!");
    // Make a copy of the data to ensure C-contiguity and prevent it going out
    // of scope
    data_ = copy_3d_array(data);
    working_=true;
    return_normals_ = return_normals;
    ready_=false;
    try {
        geom_ = std::async(std::launch::async, &Contour_Thread_Mgr::contour_surface_thread_, this, data_, threshold, cap_faces);
    } catch (...) {
        working_ = false;
        throw;
    }
}

Contour_Geometry Contour_Thread_Mgr::contour_surface_thread_(py::array_t<float> data, float threshold, bool cap_faces)
{
    auto dbuf = data.request();
    Index size[3] = { static_cast<Index>(dbuf.shape[2]),
                      static_cast<Index>(dbuf.shape[1]),
                      static_cast<Index>(dbuf.shape[0])};
    //std::cerr << "Size: " << size[0] << ", " << size[1] << ", " << size[2] << std::endl;
    int stride_size = sizeof(float);
    Stride stride[3] = { dbuf.strides[2]/stride_size, dbuf.strides[1]/stride_size, dbuf.strides[0]/stride_size };
    //std::cerr << "Strides: " << stride[0] << ", " << stride[1] << ", " << stride[2] << std::endl;
    auto cptr = std::unique_ptr<Contour_Surface>((surface((float*)dbuf.ptr, size, stride, threshold, cap_faces)));
    //std::cerr << "Vertex count: " << cptr->vertex_count() << " Triangle count: " << cptr->triangle_count() << std::endl;
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
        .def("get_result", [](Contour_Thread_Mgr& self)
        {
            auto geom = self.get_result();
            py::capsule free_vertices_when_done(geom.vertex_xyz, &delete_when_done<float>);
            py::capsule free_triangles_when_done(geom.tv_indices, &delete_when_done<int>);
            const auto& vc = geom.vertex_count;
            const auto& tc = geom.triangle_count;
            py::array_t<float> va({vc, 3}, // shape
                 {3*4, 4}, // C-style contiguous strides for n*3 float
                 geom.vertex_xyz, free_vertices_when_done );
            py::array_t<int> ta({tc, 3}, // shape
                {3*4, 4}, // C-style contiguous strides for n*3 int32
                geom.tv_indices, free_triangles_when_done );
            if (!self.return_normals())
                return py::make_tuple(va, ta);
            py::capsule free_normals_when_done(geom.normals, &delete_when_done<float>);
            py::array_t<float> na({vc, 3}, // shape
                {3*4, 4}, // C-style contiguous strides for n*3 float
                geom.normals, free_normals_when_done);
            return py::make_tuple(va, ta, na);
        })
        ;


}
