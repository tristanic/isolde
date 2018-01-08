
#include "sparse_array.h"

template <typename key_type, typename mapped_type>
sparse_array<key_type, mapped_type>::sparse_array(size_t n, key_type *keys, mapped_type *values)
{
    for (size_t i=0; i<n; ++i) {
        _data[keys[i]] = values[i];
    }
}

template <typename key_type, typename mapped_type>
sparse_array<key_type, mapped_type>::sparse_array(size_t n, key_type *keys, mapped_type *values, const mapped_type &default_value)
{
    for (size_t i=0; i<n; ++i) {
        _data[keys[i]] = values[i];
    }
    _default_value = default_value;
    _default_set = true;
}

template <typename key_type, typename mapped_type>
sparse_array<key_type, mapped_type>::sparse_array(const std::vector<key_type> &keys, const std::vector<mapped_type> &values)
{
    size_t i = 0;
    for (auto k: keys) {
        _data[k] = values[i++];
    }
}

template <typename key_type, typename mapped_type>
sparse_array<key_type, mapped_type>::sparse_array(const std::vector<key_type> &keys, const std::vector<mapped_type> &values, const mapped_type &default_value)
{
    size_t i = 0;
    for (auto k: keys) {
        _data[k] = values[i++];
    }
    _default_set = true;
    _default_value = default_value;
}

template <typename key_type, typename mapped_type>
void
sparse_array<key_type, mapped_type>::set_default_value(const mapped_type &val)
{
    _default_value = val;
    _default_set = true;
}

template <typename key_type, typename mapped_type>
const mapped_type&
sparse_array<key_type, mapped_type>::get_default_value() const
{
    if (!_default_set)
        throw std::out_of_range(err_msg_get_default_none_set());
    return _default_value;
}

template <typename key_type, typename mapped_type>
const mapped_type& 
sparse_array<key_type, mapped_type>::at(const key_type &__k) const
{
    try {
        return _data.at(__k);
    } catch (std::out_of_range) {
        if (!_default_set)
            throw std::out_of_range(err_msg_out_of_range_no_default_set());
        return _default_value;
    }
}

template <typename key_type, typename mapped_type>
const mapped_type& 
sparse_array<key_type, mapped_type>::operator[](const key_type &__k) const
{
    return this->at(__k);
}


template class sparse_array<size_t, double>;
template class sparse_array<size_t, float>;
