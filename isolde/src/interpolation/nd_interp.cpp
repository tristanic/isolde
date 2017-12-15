
#include "nd_interp.h"
#include <time.h>
#include <Python.h>

#ifdef _WIN32
# define EXPORT __declspec(dllexport)
#else
# define EXPORT __attribute__((__visibility__("default")))
#endif

RegularGridInterpolator::RegularGridInterpolator(const size_t& dim, 
        size_t* n, double* min, double* max, double* data)
{
    _dim = dim;
    size_t this_n, this_min, this_max, d_count=1;
    double step, dval;
    for (size_t i=0; i<dim; ++i) {
        this_n = n[i];
        this_min = min[i];
        this_max = max[i];
        _n.push_back(this_n);
        _min.push_back(this_min);
        _max.push_back(this_max);
        step = (this_max-this_min)/(double)(this_n-1);
        _step.push_back(step);
        _jump.push_back((size_t)pow(2,i));
        dval = this_min;
        std::vector<double> axis;
        
        for (;dval<=this_max;) {
            axis.push_back(dval);
            dval+=step;
        }
        _axes.push_back(axis);
        d_count *= this_n;
    }
    
    for (size_t i=0; i<d_count; ++i) {
        _data.push_back(data[i]);
    }
    _n_corners = (size_t)pow(2.0, (double)dim);
    corner_offsets();
} //RegularGridInterpolator


void
RegularGridInterpolator::lb_index_and_offsets(double *axis_vals, size_t &lb_index, 
    std::vector<std::pair<double, double> > &offsets)
{
//    for (size_t axis=0; axis<_dim; ++axis) {
    size_t axis_prod = 1;
    for (int axis=_dim-1; axis>=0; --axis) {
        const double &max = _max[axis];
        const double &min = _min[axis];
        const double &value = axis_vals[axis];
        if (value <= min || value >= max) {
            throw std::range_error("Value outside of interpolation range!");
        }
        size_t li = (size_t)((value-min)/_step[axis]);
        lb_index +=axis_prod*li;
        axis_prod*=_n[axis];
        const double &low = _axes[axis][li++];
        const double &high = _axes[axis][li];
        double offset = (value-low)/(high-low);
        offsets[axis]=(std::pair<double, double> (offset, 1-offset));
    }
}


/*
 * This comes out looking a little like black magic, so requires a bit
 * of explanation. We want to get the values at all the corners
 * in a well-defined order. Using the 3D case as an example, if our 
 * lower bound is (0,0,0), we want the corners in the order: 
 * ((0,0,0),(0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,0,1),(1,1,0),(1,1,1))
 * ... which is 0 to 7 in binary. The logic below simply extends this 
 * to n dimensions.
 */

void 
RegularGridInterpolator::corner_offsets()
{
    for (size_t i=0; i < _n_corners; ++i) {
        size_t corner = 0;
        size_t dim_prod = 1;
        for (size_t j=0; j<_dim; ++j) {
            corner += dim_prod * ((i & (1<<j))>>j);
            dim_prod *= _n[j];
        }
        _corner_offsets.push_back(corner);
    }
        
}


void
RegularGridInterpolator::corner_values(const size_t &lb_index, std::vector<double> &corners)
{
    for (size_t i=0; i<_corner_offsets.size(); i++) {
        corners[i]=(_data[lb_index + _corner_offsets[i]]);
    }
}

// Reduces the vector of corners in-place for efficiency
void
RegularGridInterpolator::_interpolate(const size_t &dim, std::vector<double> &corners, size_t size,
    const std::vector<std::pair<double, double> > &offsets, double* value)
{
    for (size_t i=0; i<dim; ++i) {
        const std::pair<double, double> &this_offset=offsets[dim-i-1];
        for (size_t ind=0, j=0; j<size; ind++, j+=2) {
            _interpolate1d(this_offset, corners[j], corners[j+1], &corners[ind]);
        }
        size/=2;
    }
    *value=corners[0];
}


void
RegularGridInterpolator::_interpolate1d(const std::pair<double, double> &offset, const double &lower, const double &upper, double *val)
{
    *val= offset.first*upper + offset.second*lower;
}

void RegularGridInterpolator::interpolate (double* axis_vals, const size_t &n, double* values)
{
    std::vector<std::pair<double, double>> offsets(_dim);
    std::vector<double> corners(_n_corners);
    for (size_t i=0; i<n; ++i) {
        size_t lb_index = 0;
        // find the minimum corner of the hypercube, and the offsets 
        // along each axis
        lb_index_and_offsets(axis_vals+i*_dim, lb_index, offsets);
        // ... and get values at all the corners surrounding the target
        // position.
        corner_values(lb_index, corners);
        _interpolate(_dim, corners, corners.size(), offsets, values++);

    }
    
    
    
}

//--------------------------------------------------------
// RegularGridInterpolator

// Argument delcaration types:
//
// numpy array arguments are sized, so use uint8_t for numpy's uint8,
// float32_t for numpys float32_t, etc.  The integer _t types are from
// <stdint.h>.  Special case is for numpy/C/C++ booleans which are
// processed in all cases as bytes:
//      1 == numpy.bool_().nbytes in Python
//      1 == sizeof (bool) in C++ and in C from <stdbool.h>
//      25 == sizeof (bool [25]) in C++ and C
//
// Other arguments are their normal C types and are specified with the
// appropriate ctypes annotations on the Python side.
//
// There should be very few 'int' specifications.  Any int-array should
// have a specific size, eg., int32_t, for its elements.
//
typedef uint8_t npy_bool;
typedef float float32_t;
typedef double float64_t;
typedef void *pyobject_t;

inline PyObject* unicode_from_string(const char *data, size_t size)
{
    return PyUnicode_DecodeUTF8(data, size, "replace");
}

inline PyObject* unicode_from_string(const std::string& str)
{
    return PyUnicode_DecodeUTF8(str.data(), str.size(), "replace");
}

//~ template <int len, char... description_chars>
//~ inline PyObject* unicode_from_string(const chutil::CString<len, description_chars...>& cstr)
//~ {
    //~ return PyUnicode_DecodeUTF8(static_cast<const char*>(cstr), cstr.size(),
                            //~ "replace");
//~ }

//~ inline PyObject* unicode_from_character(char c)
//~ {
    //~ char buffer[2];
    //~ buffer[0] = c;
    //~ buffer[1] = '\0';
    //~ return unicode_from_string(buffer, 1);
//~ }

static void
molc_error()
{
    // generic exception handler
    if (PyErr_Occurred())
        return;   // nothing to do, already set
    try {
        throw;    // rethrow exception to look at it
    } catch (std::bad_alloc&) {
        PyErr_SetString(PyExc_MemoryError, "not enough memory");
    } catch (std::invalid_argument& e) {
        PyErr_SetString(PyExc_TypeError, e.what());
    } catch (std::length_error& e) {
        PyErr_SetString(PyExc_MemoryError, e.what());
    } catch (std::out_of_range& e) {
        PyErr_SetString(PyExc_IndexError, e.what());
    } catch (std::overflow_error& e) {
        PyErr_SetString(PyExc_OverflowError, e.what());
    } catch (std::range_error& e) {
        PyErr_SetString(PyExc_IndexError, e.what());
    } catch (std::underflow_error& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
    } catch (std::logic_error& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
    } catch (std::ios_base::failure& e) {
        PyErr_SetString(PyExc_IOError, e.what());
    } catch (std::runtime_error& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
    } catch (std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
    } catch (...) {
        PyErr_SetString(PyExc_RuntimeError, "unknown C++ exception");
    }
}

// wrap an arbitrary function
template <typename F, typename... Args> auto
error_wrap(const F& func, Args... args) -> decltype(func(args...))
{
    try {
        return func(args...);
    } catch (...) {
        molc_error();
        return decltype(func(args...))();
    }
}

// wrap a member function
template <typename R, typename T, typename... Args> R
error_wrap(T* inst, R (T::*pm)(Args...), Args... args)
{
    try {
        return (inst->*pm)(args...);
    } catch (...) {
        molc_error();
        return R();
    }
}

// wrap a constant member function
template <typename R, typename T, typename... Args> R
error_wrap(T* inst, R (T::*pm)(Args...) const, Args... args)
{
    try {
        return (inst->*pm)(args...);
    } catch (...) {
        molc_error();
        return R();
    }
}

// wrap getting array elements via const member function
template <typename T, typename Elem, typename Elem2 = Elem> void
error_wrap_array_get(T** instances, size_t n, Elem (T::*pm)() const, Elem2* args)
{
    try {
        for (size_t i = 0; i < n; ++i)
            args[i] = (instances[i]->*pm)();
    } catch (...) {
        molc_error();
    }
}

// wrap setting array elements via member function
template <typename T, typename Elem, typename Elem2 = Elem> void
error_wrap_array_set(T** instances, size_t n, void (T::*pm)(Elem), Elem2* args)
{
    try {
        for (size_t i = 0; i < n; ++i)
            (instances[i]->*pm)(args[i]);
    } catch (...) {
        molc_error();
    }
}


extern "C" EXPORT void*
rg_interp_new(size_t dim, size_t* n, double* min, double* max, double* data)
{
    try {
        return new RegularGridInterpolator(dim, n, min, max, data);
    } catch (...) {
        molc_error();
    } return nullptr;
}

extern "C" EXPORT void
rg_interp_delete(void *ptr)
{
    RegularGridInterpolator *rg = static_cast<RegularGridInterpolator *>(ptr);    
    delete rg;
}

extern "C" EXPORT void
rg_interpolate(void* ptr, double* axis_vals, size_t n, double* values) 
{
    RegularGridInterpolator *rg = static_cast<RegularGridInterpolator *>(ptr);
    rg->interpolate(axis_vals, n, values);
}




//~ int 
//~ main() {
    //~ size_t dim = 3;
    //~ size_t n[3] {2,2,2};
    //~ double min[3] {0,0,0};
    //~ double max[3] {1,1,1};
    //~ double data[8] {0,0,0,0,0,0,1};
    //~ RegularGridInterpolator interp(dim, n, min, max, data);
    //~ double* results = new double[100000];
    //~ //double axis_vals[9] {0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.8, 0.8, 0.8};
    //~ double* axis_vals = new double[300000];
    //~ for (size_t i=0; i<300000; ++i)
        //~ axis_vals[i]=(double)rand()/RAND_MAX*0.95+0.01;
    
    //~ double start = (double)clock();
    
    //~ //for (size_t i=0; i<100000; ++i) {
        //~ interp.interpolate(axis_vals, 100000, results);
    //~ //}
    
    //~ std::cout << "100,000 3D interpolations took " << ((double)clock()-start)/CLOCKS_PER_SEC << " seconds." << std::endl;
    //~ for (size_t i=0; i<5; ++i) {
        //~ std::cout << results[i] << std::endl;
    //~ }
    //~ delete results;
    //~ delete axis_vals;
    
    //~ return 0;
//~ }
