#ifndef ISOLDE_COLOR
#define ISOLDE_COLOR

#include <cstdint>
#include <algorithm>

namespace isolde
{

namespace colors
{

//! RGBA colours in range 0..1
typedef double color[4];

//! RGBA colours in range 0..255
typedef uint8_t intcolor[4];

inline void intcolor_as_color(intcolor &in, color &out)
{
    for (size_t i=0; i<4; ++i)
        out[i] = (double)in[i] / 255.0;
}

inline void color_as_intcolor(color &in, intcolor &out)
{
    for (size_t i=0; i<4; ++i)
        out[i] = (uint8_t)(in[i] * 255.0);
}


static const color grey = {0.5, 0.5, 0.5, 1.0};

inline void copy_color(color &from, color &to)
{
    for (size_t i=0; i<4; ++i) {
        to[i] = from[i];
    }
}

inline void interpolate_colors(const color &min_color, const color &max_color,
    const double &min_val, const double &max_val, const double &score, color &out)
{
    double offset = (score-min_val)/(max_val-min_val);
    for (size_t i=0; i<4; ++i) {
        out[i] = min_color[i] + (max_color[i]-min_color[i])*offset;
    }    
}


class colormap
{
private:
    std::unordered_map<double, color> _cmap;
    std::vector<double> _mapped_vals;
    color _above_max;
    color _below_min;
    bool _color_above = false;
    bool _color_below = false;
public:    
    colormap() {}
    ~colormap() {}
    colormap(double* mapped_vals, color* mapped_colors, size_t n)
        : _color_above(false), _color_below(false)
    {
        for (size_t i=0; i<n; ++i)
        {
            auto &mcolor = _cmap[mapped_vals[i]];
            for (size_t j=0; j<4; ++j)
                mcolor[j] = mapped_colors[i][j];
            _mapped_vals.push_back(mapped_vals[i]);
        }
        std::sort(_mapped_vals.begin(), _mapped_vals.end());
    }
    colormap(double* mapped_vals, color* mapped_colors, size_t n,
        const color &above_max, const color &below_min)
        : colormap(mapped_vals, mapped_colors, n)
    {
        for(size_t i=0; i<4; ++i)
        {
            _above_max[i] = above_max[i];
            _below_min[i] = below_min[i];
            _color_above = true;
            _color_below = true;
        }
    }
    
    void interpolate(double value, color &rgba)
    {
        interpolate(&value, 0, &rgba);
    }

    void interpolate(double *values, size_t n, color *rgba)
    {
        color this_color;
        for (size_t i=0; i<n; ++i)
        {
            double v = values[i];
            if (v<_mapped_vals[0]) {
                if (_color_below) {
                    copy_color(_below_min, this_color);
                } else {
                    copy_color((_cmap[_mapped_vals[0]]), this_color);
                }
            } else if (v>_mapped_vals.back()) {
                if (_color_above) {
                    copy_color(_above_max, this_color);
                } else {
                    copy_color(_cmap[_mapped_vals.back()], this_color);
                }
            } else {
                // first value greater than or equal to v
                auto gt = std::lower_bound(_mapped_vals.begin(), _mapped_vals.end(), v);
                auto lt = gt-1;
                interpolate_colors(_cmap[*lt], _cmap[*gt], *lt, *gt, v, this_color);
            }
            copy_color(this_color, rgba[i]);
        }
    }
};


} //namespace colors
} //namespace isolde
#endif //ISOLDE_COLOR
