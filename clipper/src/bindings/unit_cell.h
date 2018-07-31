#pragma once

#include <set>

#include <clipper/clipper.h>
#include "symops.h"
namespace clipper
{

/*
  Used for finding symmetry operator(s) mapping a set of coordinates back to a
  region in the unit cell. For internal use only; no need to be exposed to
  Python.
*/
struct Symop_and_Unit_Cell_Offset
{
  int symop;
  Coord_grid uc_offset;

  Symop_and_Unit_Cell_Offset(int op, const Coord_grid& offset)
  {
    symop = op;
    uc_offset = offset;
  }

  bool operator<(const Symop_and_Unit_Cell_Offset& other) const
  {
    if (symop < other.symop ) return true;
    if (symop > other.symop ) return false;
    const Coord_grid& uo = other.uc_offset;
      for (int i=0; i<3; ++i) {
        if (uc_offset[i] < uo[i]) return true;
        if (uc_offset[i] > uo[i]) return false;
      }
    return false;
  }

  bool operator==(const Symop_and_Unit_Cell_Offset& other) const
  {
    if (symop != other.symop) return false;
    const Coord_grid& uo = other.uc_offset;
    for (int i=0; i<3; ++i) {
      if (uc_offset[i] != uo[i]) return false;
    }
    return true;
  }

};


class Unit_Cell
{
public:
    ~Unit_Cell() {};

    Unit_Cell(const Coord_frac& ref, const Atom_list& atoms, const Cell& cell,
        const Spacegroup& sg, const Grid_sampling& grid, int padding=0);

    inline const Symops& symops() const { return symops_; } // Getter
    // inline Symops& symops() { return symops_; } // Setter

    inline const Symops& inv_symops() const { return inv_symops_; } // Getter
    // inline Symops& inv_symops() { return inv_symops(); } // Setter

    inline const Isymops& isymops() const { return isymops_; } // Getter
    // inline Isymops& isymops() const { return isymops_; } // Setter

    inline const Isymops& inv_isymops() const { return inv_isymops_; } // Getter
    // inline Isymops& inv_isymops() { return inv_isymops_; } // Setter

    inline const Coord_frac& ref_coord() const { return ref_; } // Getter
    inline void set_ref_coord(const Coord_frac& new_ref) {
        ref_ = new_ref;
        recalculate_symops_();
    }

    inline const Cell& cell() const { return cell_; }
    inline const Spacegroup& spacegroup() const { return sg_; }
    inline const Grid_sampling& grid() const { return grid_; }


    inline const Grid_range& ref_box() const { return reference_model_bounds_; }

    inline size_t ref_box_shortest_side() const { return min_side_; }

    inline const Coord_grid& min() const { return cell_origin_; }
    inline Coord_grid max() const {
        return cell_origin_ + Coord_grid(grid_.nu()-1, grid_.nv()-1, grid_.nw()-1); }

    void find_symops_for_coord(std::set<Symop_and_Unit_Cell_Offset>& pairlist,
        const Coord_grid& coord) const;

    Symops all_symops_in_box(const Grid_range& range, bool always_include_identity=false, int sample_frequency=2) const;
    Symops all_symops_in_box(const Coord_orth& origin_xyz, const Vec3<int>& box_size_uvw,
        bool always_include_identity=false, int sample_frequency=2) const;


    /* Find the minimum and maximum grid coordinates of a box encompassing the atoms,
     * pad it by padding grid units either side, and save the result as a Grid_range.
     */
    void update_reference_model_bounds(const Atom_list& atoms, const int& padding);


private:

    Cell cell_;
    Spacegroup sg_;
    Grid_sampling grid_;
    Coord_frac ref_;  // Our reference coordinate in fractional...
    Coord_grid grid_ref_; // ... and grid coordinates.
    Symops symops_; // The list of RTops mapping our reference coordinate to equivalent positions in the same unit cell
    Symops inv_symops_; // For convenience, the inverse symops are also saved.
    Isymops isymops_; // Integerised symops for fast operations on grid coordinates
    Isymops inv_isymops_; //... and their inverses
    std::vector<int> symop_map_; // Mapping of the indices of symops_ to Clipper's internal symop indices
    //Grid_range ref_asu_;
    Coord_grid cell_origin_; // The origin of the unit cell containing our reference model
    Grid_range reference_model_bounds_; // Smallest grid range encompassing our atomic model
    int min_side_; // The length of the smallest side of the grid box in reference_model_bounds
    /* Because life is never easy, the reference model will almost certainly
     * span across multiple unit cells (unlike the nice, neat map reference
     * asu used by Clipper. Therefore, to find symmetry equivalents we
     * need to test the main unit cell plus all 26 direct neighbours.
     */
    std::vector<Coord_grid> ref_model_cell_origins_;

    // If the reference coordinate is changed, we need to re-calculate the
    // symop arrays to make sure they all still map back to this unit cell.
    void recalculate_symops_();
}; // class Unit_Cell


// Implementations

Unit_Cell::Unit_Cell(const Coord_frac& ref, const Atom_list& atoms, const Cell& cell,
    const Spacegroup& sg, const Grid_sampling& grid, int padding)
    : cell_(cell),  sg_(sg), grid_(grid), ref_(ref)
{
    /* Store the nearest grid coordinate to the input reference coordinate.
     * Make the equivalent grid, fractional and map reference coords
     * and store for repeated use.
     */
    grid_ref_ = ref.coord_grid(grid);
    cell_origin_ = grid_ref_ - grid_ref_.unit(grid);
    ref_ = grid_ref_.coord_frac(grid);
    update_reference_model_bounds(atoms, padding);

    // In order to catch all corner cases, we'll need to search the 3x3 block
    // of unit cells consisting of this one and its direct neighbours.
    Coord_grid this_offset;
    for (int i = -1; i < 2; i++) {
      for (int j = -1; j < 2; j++) {
        for (int k = -1; k < 2; k++) {
          this_offset = Coord_grid(grid.nu()*i, grid.nv()*j, grid.nw()*k);
          ref_model_cell_origins_.push_back(cell_origin_ + this_offset);
        }
      }
    }
    recalculate_symops_();

}

void
Unit_Cell::recalculate_symops_()
{
    // For each symop in this spacegroup, we need to find the whole-unit-cell
    // translations needed to bring the transformed reference coordinate back
    // into this unit cell. The result will be stored as a Symops array.
    Coord_frac origin_frac = cell_origin_.coord_frac(grid_);
    symops_ = Symops();
    inv_symops_ = Symops();
    for (int i=0; i<sg_.num_symops(); ++i)
    {
        const Symop& this_symop = sg_.symop(i);
        Coord_frac tc = this_symop*ref_;
        Coord_frac cell_offset = tc - tc.lattice_copy_unit() - origin_frac;
        symops_.append(this_symop, -cell_offset);
    }

    // Pre-calculate the inverse and integer symops
    for (size_t i=0; i<symops_.size(); ++i)
    {
        RTop_frac inv_op = symops_.with_cell_translation(i).inverse();
        inv_symops_.append(inv_op);
    }

    isymops_ = Isymops(symops_, grid_);
    inv_isymops_ = Isymops(inv_symops_, grid_);
}

void
Unit_Cell::update_reference_model_bounds(const Atom_list& atoms, const int& padding)
{
    Coord_grid ref_min = atoms[0].coord_orth().coord_frac(cell_).coord_grid(grid_);
    Coord_grid ref_max = ref_min;
    Coord_grid pad(padding,padding,padding);
    for (const auto& atom: atoms) {
      const Coord_grid& thiscoord = atom.coord_orth().coord_frac(cell_).coord_grid(grid_);
      for (size_t i = 0; i < 3; i++) {
        if (thiscoord[i] < ref_min[i]) ref_min[i] = thiscoord[i];
        else if (thiscoord[i] > ref_max[i]) ref_max[i] = thiscoord[i];
      }
    }
    reference_model_bounds_ = Grid_range( ref_min-pad, ref_max+pad );

    /* Find the side of the reference box covering the shortest
    * distance. When searching a box for relevant symops, we'll
    * sample it in steps of half this size.
    */

    size_t nu, nv, nw;
    nu = reference_model_bounds_.nu();
    nv = reference_model_bounds_.nv();
    nw = reference_model_bounds_.nw();
    min_side_ = nv < nu ? (nw < nv ? nw : nv) : nu;

}

void
Unit_Cell::find_symops_for_coord(std::set<Symop_and_Unit_Cell_Offset>& pairlist,
        const Coord_grid& coord) const
{
    Coord_grid shift, t_coord;
    for (const auto& co: ref_model_cell_origins_)
    {
        // First work out how many unit cells we are away from our reference
        shift = coord - coord.unit(grid_) - co;
        /* Now try the inverse symops, and work out which (if any) put us
         * within the bounds of our reference box. NOTE: since our box
         * is the limits of the model (not an ASU), it's entirely possible
         * for the coordinate to be in more than one box - or, indeed, no box
         * at all. Therefore, we should err on the side of caution and test
         * all symops for each point.
        */
        for (size_t i=0; i< inv_symops_.size(); ++i) {
            t_coord = (coord.unit(grid_) + co).transform(inv_isymops_.with_cell_translation(i));
            if (reference_model_bounds_.in_grid(t_coord))
                pairlist.insert(Symop_and_Unit_Cell_Offset(i, shift));
        }
    }
}


Symops
Unit_Cell::all_symops_in_box(const Coord_orth& origin_xyz, const Vec3<int>& box_size_uvw,
    bool always_include_identity, int sample_frequency) const
{
    Coord_grid grid_min = origin_xyz.coord_frac(cell_).coord_grid(grid_);
    Coord_grid grid_max = grid_min + Coord_grid(box_size_uvw);
    Grid_range range(grid_min, grid_max);
    return all_symops_in_box(range, always_include_identity, sample_frequency);
}

Symops
Unit_Cell::all_symops_in_box(const Grid_range& range, bool always_include_identity, int sample_frequency) const
{
    static const Symop_and_Unit_Cell_Offset identity(0, Coord_grid(0,0,0));
    Symops ret;
    std::set<Symop_and_Unit_Cell_Offset> ops_and_translations;

    const Coord_grid& box_origin = range.min();
    const Coord_grid& box_max = range.max();

    /*
        Sample the box in steps equal to 1/frequency times the length of the
        shortest side of the box encompassing the atomic model, making sure we
        also capture the faces and corners.
    */
    int step_size = std::max(min_side_/sample_frequency, 1);
    bool u_done = false, v_done = false, w_done = false;
    int u,v,w;
    Coord_grid thiscg;
    u = box_origin[0];
    while (!u_done) {
        v=box_origin[1];
        v_done=false;
        if (u == box_max[0]) u_done=true;
        while (!v_done) {
            w = box_origin[2];
            w_done = false;
            if (v==box_max[1]) v_done = true;
            while (!w_done) {
                if (w == box_max[2]) w_done = true;
                thiscg.u()=u; thiscg.v()=v; thiscg.w()=w;
                find_symops_for_coord(ops_and_translations, thiscg);
                w = std::min(w+step_size, box_max[2]);
            }
            v = std::min(v+step_size, box_max[1]);
        }
        u = std::min(u+step_size, box_max[0]);
    }

    auto it = ops_and_translations.find(identity);
    if (it != ops_and_translations.end()) {
        ops_and_translations.erase(it);
        ret.append(RTop_frac::identity());
    } else if (always_include_identity)
        ret.append(RTop_frac::identity());

    for (const auto& ot: ops_and_translations) {
        RTop_frac thisop = symops_.with_cell_translation(ot.symop);
        thisop.trn() += ot.uc_offset.coord_frac(grid_);
        ret.append(thisop);
    }
    return ret;
}


} // namespace clipper
