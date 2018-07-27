#pragma once
#include <clipper/clipper.h>

namespace clipper
{

// Get the number of whole cell steps in each direction for a RTop_frac
inline Coord_frac cell_shift (const Coord_frac& coord)
{
    Coord_frac ret;
    ftype intpart;
    for (size_t i = 0; i < 3; i++) {
        modf(coord[i], &intpart);
        ret[i] = intpart;
    }
    return ret;
}

inline Coord_grid cell_shift (const Coord_grid& coord, const Grid_sampling& grid)
{
    Coord_grid ret;
    ret = coord - coord.unit(grid);
    return ret;
}

class Symops
{
public:
    Symops() {} // Default constructor
    inline Symops(const std::vector<RTop_frac>& ops)
    {
        for (const auto& op: ops) {
            append(op);
        }
    }
    inline Symops(const std::vector<Symop>& ops)
    {
        for (const auto& op: ops) {
            append(op, Coord_frac(0,0,0));
        }
    }

    inline void append(const Symop& op, const Coord_frac& offset)
    {
        RTop_frac with_trn = RTop_frac(op);
        with_trn.trn() += static_cast<Vec3<ftype>>(offset);
        symops_.push_back(op);
        unit_translations_.push_back(offset);
        symops_with_cell_translations_.push_back(with_trn);
    }

    inline void append(const RTop_frac& op)
    {
        Symop symop (op);
        auto unit_trn = cell_shift(Coord_frac(op.trn()));
        RTop_frac with_trn = RTop_frac(symop);
        with_trn.trn() += static_cast<Vec3<ftype>>(unit_trn);
        symops_.push_back(symop);
        unit_translations_.push_back(unit_trn);
        symops_with_cell_translations_.push_back(with_trn);
    }

    // getter
    inline const Symop& at(const int& i) const { return symops_.at(i); }

    // setter
    inline void replace(const int& i, const RTop_frac& op)
    {
        Symop symop(op);
        auto unit_trn = cell_shift(Coord_frac(op.trn()));
        RTop_frac with_trn = RTop_frac(symop);
        with_trn.trn() += static_cast<Vec3<ftype>>(unit_trn);
        symops_.at(i) = symop;
        unit_translations_[i] = unit_trn;
        symops_with_cell_translations_[i] = with_trn;
    }

    //! get element without bounds checking
    inline const Symop& operator []( const int& i ) const { return symops_[i]; }

    inline RTop_frac pop(int i)
    {
        RTop_frac ret = symops_.at(i);
        symops_.erase(symops_.begin()+i);
        unit_translations_.erase(unit_translations_.begin()+i);
        symops_with_cell_translations_.erase(symops_with_cell_translations_.begin()+i);
        return ret;
    }

    inline size_t size() const { return symops_.size(); }

    // Return the integer component of the translation for each RTop_frac.
    // For use in C++ layer only (no bounds checking).
    inline const Coord_frac& cell_trans( const int& i) const
    {
        return unit_translations_[i];
    }

    inline const RTop_frac& with_cell_translation(const int& i) const
    {
        return symops_with_cell_translations_.at(i);
    }

private:
    std::vector<Symop> symops_;
    std::vector<Coord_frac> unit_translations_;
    std::vector<RTop_frac> symops_with_cell_translations_;
}; // class Symops

// Not exposed to python.
class Isymops
{
public:
    Isymops() {}
    inline Isymops(const Symops& ops, const Grid_sampling& grid)
    {
        grid_ = grid;
        for (size_t i=0; i<ops.size(); ++i) {
            Isymop op = Isymop(ops[i], grid);
            isymops_.push_back(op);
            Coord_grid trn = ops.cell_trans(i).coord_grid(grid);
            unit_translations_.push_back(trn);
            op.trn()+= trn;
            isymops_with_cell_translations_.push_back(op);
        }
    }
    ~Isymops() {}

    //! get element without bounds checking
    inline const Isymop& operator []( const int& i ) const { return isymops_[i]; }

    inline const Isymop& at(const size_t& i) const { return isymops_.at(i); }
    inline size_t size() const { return isymops_.size(); }

    inline void append(const Isymop& op, const Coord_grid& cell_trn)
    {
        isymops_.push_back(op);
        unit_translations_.push_back(cell_trn);
        Isymop with_trn(op);
        with_trn.trn() += static_cast<Vec3<int>>(cell_trn);
        isymops_with_cell_translations_.push_back(with_trn);
    }

    inline Isymop pop(const int& i)
    {
        auto ret = isymops_.at(i);
        isymops_.erase(isymops_.begin()+i);
        unit_translations_.erase(unit_translations_.begin()+i);
        return ret;
    }

    // No bounds checking.
    inline const Coord_grid& cell_trans(const int& i) const
    {
        return unit_translations_[i];
    }

    inline const Isymop& with_cell_translation(const int& i) const
    {
        return isymops_with_cell_translations_.at(i);
    }


private:
    std::vector<Isymop> isymops_;
    std::vector<Coord_grid> unit_translations_;
    std::vector<Isymop> isymops_with_cell_translations_;
    Grid_sampling grid_;
}; // class Isymops



} // namespace clipper
