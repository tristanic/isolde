
#include "nd_interp.h"
#include <time.h>


#ifdef _WIN32
# define EXPORT __declspec(dllexport)
#else
# define EXPORT __attribute__((__visibility__("default")))
#endif

RegularGridInterpolator::RegularGridInterpolator(const size_t& dim, 
        size_t* n, double* min, double* max, double* data)
{
    _dim = dim;
    size_t this_n, d_count=1;
    double this_min, this_max;
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

extern "C"
{

EXPORT void*
rg_interp_new(size_t dim, size_t* n, double* min, double* max, double* data)
{
    try {
        return new RegularGridInterpolator(dim, n, min, max, data);
    } catch (...) {
        molc_error();
    } return nullptr;
}

EXPORT void*
rg_interp_copy(void *ptr)
{
    RegularGridInterpolator *rg = static_cast<RegularGridInterpolator *>(ptr);
    return new RegularGridInterpolator(*rg);
}

EXPORT void
rg_interp_delete(void *ptr)
{
    try {
        RegularGridInterpolator *rg = static_cast<RegularGridInterpolator *>(ptr);
        delete rg;
    } catch (...) {
        molc_error();
        return;
    }   
}

EXPORT void
rg_interpolate(void* ptr, double* axis_vals, size_t n, double* values) 
{
    try {
        RegularGridInterpolator *rg = static_cast<RegularGridInterpolator *>(ptr);
        rg->interpolate(axis_vals, n, values);
    } catch (...) {
        molc_error();
        return;
    }
}

EXPORT void
rg_interp_min(void* ptr, double* ret)
{
    try {
        RegularGridInterpolator *rg = static_cast<RegularGridInterpolator *>(ptr);
        for (auto m: rg->min()) {
            *ret++ = m;
        }
    } catch (...) {
        molc_error();
        return;
    }
}

EXPORT void
rg_interp_max(void* ptr, double* ret)
{
    try {
        RegularGridInterpolator *rg = static_cast<RegularGridInterpolator *>(ptr);
        for (auto m: rg->max()) {
            *ret++ = m;
        }
    } catch (...) {
        molc_error();
        return;
    }
}

EXPORT void
rg_interp_lengths(void* ptr, size_t* ret)
{
    try {
        RegularGridInterpolator *rg = static_cast<RegularGridInterpolator *>(ptr);
        for (auto l: rg->length()) {
            *ret++ = l;
        }
    } catch (...) {
        molc_error();
        return;
    }
}



EXPORT size_t
rg_interp_dim(void* ptr)
{
    try {
        RegularGridInterpolator *rg = static_cast<RegularGridInterpolator *>(ptr);
        return rg->dim();
    } catch (...) {
        molc_error();
        return 0;
    }
}

EXPORT void
rg_interp_values(void* ptr, double* ret)
{
    try {
        RegularGridInterpolator *rg = static_cast<RegularGridInterpolator *>(ptr);
        for (auto d: rg->data()) {
            *ret++ = d;
        }
    } catch (...) {
        molc_error();
        return;
    }
}




} // extern "C"


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
