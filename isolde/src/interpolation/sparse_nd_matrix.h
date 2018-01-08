
#ifndef SPARSE_ARRAY
#define SPARSE_ARRAY

#include <unordered_map>
#include <stdint.h>
#include <vector>
#include <stdexcept>


/**
 * @brief A sparse array implementation
 * 
 * Uses a std::unordered_map to emulate a sparse 1D array, where 
 * accessing elements that have not yet been assigned returns a 
 * default value. Intended for interpolation on high-dimensional maps
 * where most of the values are zero.
 */
template <typename key_type, typename mapped_type>
class sparse_array
{
public:
    sparse_array(size_t n, key_type *keys, mapped_type *values);
    sparse_array(size_t n, key_type *keys, mapped_type *values, const mapped_type &default_value);
    sparse_array(const std::vector<key_type> &keys, const std::vector<mapped_type> &values);
    sparse_array(const std::vector<key_type> &keys, const std::vector<mapped_type> &values, const mapped_type &default_value);
    ~sparse_array();
    mapped_type &default() {return _default_value;} //setter
    const mapped_type &default() const {return _default_value;} //getter
    /**
     * @brief Setter. Simply wraps the std::unordered_map [] operator.
     * 
     * DO NOT USE THIS FOR ACCESSING VALUES. Using an index not currently
     * in the array will add some (possibly undefined) default mapped_type
     * entry. Use at() instead.
     */
    mapped_type& operator[](const key_type &__k) {return _data[__k];};
    const mapped_type& at(const key_type &__k) const;
    

private:
    std::unordered_map<key_type, mapped_type> _data;
    mapped_type _default_value;
    
    
    
}; //class sparse_array

#endif //SPARSE_ARRAY
