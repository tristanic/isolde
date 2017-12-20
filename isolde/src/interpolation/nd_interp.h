//! N-dimensional regular grid interpolator

#ifndef ND_GRID_INTERP
#define ND_GRID_INTERP

#include <stdint.h>
#include <vector>
#include <math.h>
#include <stdexcept>
#include <iostream>
#include "../molc.h"

template <typename T>
class RegularGridInterpolator
{

public:
    RegularGridInterpolator(); // null constructor
    ~RegularGridInterpolator() {} // destructor
    //! Construct a RegularGridInterpolator object for the given data
    /*!
     * This implementation requires one data value for every grid point
     * dim: the number of dimensions
     * n:   the number of posize_ts in each dimension
     * min: the minimum axis value for each dimension
     * max: the maximum axis value for each dimension
     * data: the actual data to be interpolated (must match the dimensions
     *       defined by the previous arguments)
     */
    RegularGridInterpolator(const size_t &dim, size_t* n, T* min, T* max, T* data);
    
    
    //! Construct a RegularGridInterpolator object for the given data
    /*!
     * This implementation creates a sparse n-dimensional matrix 
     * to reduce memory usage for large arrays with lots of zero values.
     * dim: the number of dimensions
     * n:   the number of posize_ts in each dimension
     * min: the minimum axis value for each dimension
     * max: the maximum axis value for each dimension
     * data_coords: the (x1,x2,x3,...,xn) coordinates for each point in data
     * data: the actual data to be interpolated (must match the number
     *       of coordinates in data_coords)
     */
    RegularGridInterpolator(const size_t &dim, size_t* n, T* min, T* max, T* data_coords, T* data);
    
    
    void interpolate(T* axis_vals, const size_t &n, T* values);
    const std::vector<T> &min() const {return _min;}
    const std::vector<T> &max() const {return _max;}
    const size_t &dim() const {return _dim;}
    const std::vector<T> &data() const {return _data;}
    const std::vector<size_t> &length() const {return _n;} 
    
    
private:
    void corner_values(const size_t &lb_indices, std::vector<T> &corners);
    void lb_index_and_offsets(T *axis_vals, size_t &lb_index, 
        std::vector<std::pair<T, T> > &offsets);
    void _interpolate(const size_t &dim, std::vector<T> &corners, size_t size,
    const std::vector<std::pair<T, T> > &offsets, T *value);
    void _interpolate1d(const std::pair<T, T> &offset, const T& lower, const T& upper, T *val);
    void corner_offsets();

    size_t _dim;
    size_t _n_corners;
    std::vector<size_t> _n;
    std::vector<T> _min;
    std::vector<T> _max;
    std::vector<T> _step;
    std::vector<std::vector<T> > _axes;
    
    //TODO: Replace _data with a std::unordered_map sparse array implementation
    //      to minimise memory use for higher dimensions
    std::vector<T> _data;
    std::vector<size_t> _corner_offsets;
    std::vector<size_t> _jump;
    
}; //RegularGridInterpolator


#endif // ND_GRID_INTERP
