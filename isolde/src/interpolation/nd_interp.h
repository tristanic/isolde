//! N-dimensional regular grid interpolator

#ifndef ND_GRID_INTERP
#define ND_GRID_INTERP

#include <stdint.h>
#include <vector>
#include <math.h>
#include <stdexcept>
#include <iostream>
#include "../molc.h"

class RegularGridInterpolator
{

public:
    RegularGridInterpolator(); // null constructor
    ~RegularGridInterpolator() {} // destructor
    //! Construct a RegularGridInterpolator object for the given data
    /*!
     * dim: the number of dimensions
     * n:   the number of posize_ts in each dimension
     * min: the minimum axis value for each dimension
     * max: the maximum axis value for each dimension
     * data: the actual data to be interpolated (must match the dimensions
     *       defined by the previous arguments)
     */
    RegularGridInterpolator(const size_t& dim, size_t* n, double* min, double* max, double* data);
    void interpolate(double* axis_vals, const size_t &n, double* values);
    const std::vector<double> &min() const {return _min;}
    const std::vector<double> &max() const {return _max;}
    const size_t &dim() const {return _dim;}
    const std::vector<double> &data() const {return _data;}
    const std::vector<size_t> &length() const {return _n;} 
    
    
private:
    void corner_values(const size_t &lb_indices, std::vector<double> &corners);
    void lb_index_and_offsets(double *axis_vals, size_t &lb_index, 
        std::vector<std::pair<double, double> > &offsets);
    void _interpolate(const size_t &dim, std::vector<double> &corners, size_t size,
    const std::vector<std::pair<double, double> > &offsets, double *value);
    void _interpolate1d(const std::pair<double, double> &offset, const double& lower, const double& upper, double *val);
    void corner_offsets();

    size_t _dim;
    size_t _n_corners;
    std::vector<size_t> _n;
    std::vector<double> _min;
    std::vector<double> _max;
    std::vector<double> _step;
    std::vector<std::vector<double> > _axes;
    std::vector<double> _data;
    std::vector<size_t> _corner_offsets;
    std::vector<size_t> _jump;
    
}; //RegularGridInterpolator


#endif // ND_GRID_INTERP
