
#ifndef SPARSE_ARRAY
#define SPARSE_ARRAY

#include <unordered_map>
#include <stdint.h>
#include <vector>
#include <stdexcept>
#include <iostream>

#ifdef _WIN32
# define EXPORT __declspec(dllexport)
#else
# define EXPORT __attribute__((__visibility__("default")))
#endif



/**
 * @brief A sparse array implementation
 * 
 * Uses a std::unordered_map to emulate a sparse 1D array, where 
 * accessing elements that have not yet been assigned returns a 
 * default value. Intended for interpolation on high-dimensional maps
 * where most of the values are zero.
 */
template <typename key_type, typename mapped_type>
class EXPORT
sparse_array
{
public:

    typedef std::unordered_map<key_type, mapped_type> _Umap;
    /**
     * Constructor 
     */
    sparse_array() {};
    sparse_array(size_t n, key_type *keys, mapped_type *values);
    sparse_array(size_t n, key_type *keys, mapped_type *values, const mapped_type &default_value);
    sparse_array(const std::vector<key_type> &keys, const std::vector<mapped_type> &values);
    sparse_array(const std::vector<key_type> &keys, const std::vector<mapped_type> &values, const mapped_type &default_value);
    ~sparse_array() {}
    void set_default_value(const mapped_type &val); //setter
    const mapped_type& get_default_value() const; //getter
    void set(const key_type& __k, const mapped_type& m) { 
        std::cerr << "Setting " << __k << " to " << m << std::endl;
        _data[__k] = m; }
    const mapped_type& operator[](const key_type &__k) const;
    const mapped_type& at(const key_type &__k) const;

private:
    const char* err_msg_get_default_none_set() const { return "No default value set!"; }
    const char* err_msg_out_of_range_no_default_set() const
    {
        return "This entry has no value, and no default has been set!";
    }
    
    _Umap _data;
    bool _default_set = false;
    mapped_type _default_value;
    
    
    
}; //class sparse_array






#endif //SPARSE_ARRAY
