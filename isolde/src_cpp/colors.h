/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 02-Apr-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#ifndef ISOLDE_COLOR
#define ISOLDE_COLOR

#include <cstdint>
#include <algorithm>
#include <vector>
#include <array>
#include <iostream>

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

template <typename T>
inline void copy_color(const T &from, color &to)
{
    for (size_t i=0; i<4; ++i) {
        to[i] = from[i];
    }
}

inline void copy_color(const intcolor &from, intcolor &to)
{
    for (size_t i=0; i<4; ++i) {
        to[i] = from[i];
    }
}

template <typename T>
inline void interpolate_colors(const T &min_color, const T &max_color,
    const double &min_val, const double &max_val, const double &score, color &out)
{
    double offset = (score-min_val)/(max_val-min_val);
    for (size_t i=0; i<4; ++i) {
        out[i] = min_color[i] + (max_color[i]-min_color[i])*offset;
    }
}

/*! Unlike colormap, variable_colormap requires the cutoff values to be provided
 *  with each call, and does not provide outside-of-range colouring.
 */
class variable_colormap
{
public:
    variable_colormap() {}
    ~variable_colormap() {}
    variable_colormap(color* mapped_colors, size_t n)
    {
        for (size_t i=0; i<n; ++i) {
            std::array<double, 4> thiscolor;
            std::copy(mapped_colors[i], mapped_colors[i]+4, thiscolor.begin());
            _colors.push_back(thiscolor);
        }
        _num_colors = n;
    }
    const std::vector<std::array<double, 4>>& mapped_colors() const { return _colors; }

    void interpolate(double value, const double *cutoffs, color &rgba)
    {
        if (value <= cutoffs[0])
            copy_color(_colors[0], rgba);
        else if (value >= cutoffs[_num_colors-1]) {
            copy_color(_colors.back(), rgba);
        }
        else {
            size_t i=0;
            for (;i<_num_colors;) {
                if ((i<_num_colors-1) && cutoffs[i] == cutoffs[i+1]) {
                    // avoid divide-by-zero
                    i++;
                    continue;
                }
                if (value < cutoffs[i])
                    break;
                i++;
            }
            interpolate_colors(_colors[i-1], _colors[i], cutoffs[i-1], cutoffs[i], value, rgba);
        }
    }
private:
    std::vector<std::array<double, 4>> _colors;
    size_t _num_colors;
};


class colormap
{
private:
    // std::unordered_map<double, color> _cmap;
    // std::vector<double> _mapped_vals;
    struct mapped_color
    {
        double val;
        color thecolor;
        bool operator< (const mapped_color &rhs) const {return val < rhs.val;}
        bool operator< (const double &v) const {return val < v;}
        bool operator> (const double &v) const {return val > v;}
        bool operator<= (const double &v) const {return val <= v;}
        bool operator>= (const double &v) const {return val >= v;}
        bool operator== (const double &v) const {return val == v;}
        bool operator!= (const double &v) const {return val != v;}



        mapped_color() {}
        mapped_color(const double &v, const color &c): val(v)
        {
            for (size_t i=0; i<4; ++i) {
                thecolor[i] = c[i];
            }
        }
    };
    std::vector<mapped_color> _mapped_colors;
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
            _mapped_colors.emplace_back(*mapped_vals++, *mapped_colors++);
        }
        std::sort(_mapped_colors.begin(), _mapped_colors.end());
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

    const std::vector<mapped_color> &mapped_colors() {return _mapped_colors;}
    const color &below_min_color() {return _below_min;}
    const color &above_max_color() {return _above_max;}

    void interpolate(double value, color &rgba)
    {
        interpolate(&value, 1, &rgba);
    }

    void interpolate(double *values, size_t n, color *rgba)
    {
        color this_color;
        for (size_t i=0; i<n; ++i)
        {
            double v = values[i];
            if (_mapped_colors[0]>v) {
                if (_color_below) {
                    copy_color(_below_min, this_color);
                } else {
                    copy_color((_mapped_colors[0].thecolor), this_color);
                }
            } else if (_mapped_colors.back()<v) {
                if (_color_above) {
                    copy_color(_above_max, this_color);
                } else {
                    copy_color(_mapped_colors.back().thecolor, this_color);
                }
            } else {
                // first value greater than or equal to v
                auto gt = std::lower_bound(_mapped_colors.begin(), _mapped_colors.end(), v);
                auto lt = gt-1;
                interpolate_colors((*lt).thecolor, (*gt).thecolor, (*lt).val, (*gt).val, v, this_color);
            }
            copy_color(this_color, rgba[i]);
        }
    }
};


} //namespace colors
} //namespace isolde
#endif //ISOLDE_COLOR
