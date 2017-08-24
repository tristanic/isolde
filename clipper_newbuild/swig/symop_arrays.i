%{
  namespace clipper
  {
  /* Some simple functions for sorting and removing duplicates from
   * vectors of (symop index, Coord_frac) or (symop index, Coord_grid) pairs.
   * These are only used internally for finding symmetry operators mapping
   * the atomic model to a region, and don't need to be wrapped.
   */
  bool compare_int_Coord_frac_pairs (const std::pair<int, Coord_frac > &a, const std::pair<int, Coord_frac > &b)
  {
    if (a.first < b.first) return true;
    if (a.first > b.first) return false;
    for (int i = 0; i < 3; i++) {
      if (a.second[i] > b.second[i]) return false;
      if (a.second[i] < b.second[i]) return true;
    }
    return false;
  }

  bool int_Coord_frac_pairs_equal (const std::pair<int, Coord_frac > &a, const std::pair<int, Coord_frac > &b)
  {
    double eps = 1e-6;
    if (a.first != b.first) return false;
    for (int i=0; i < 3; i++) {
      if (std::abs(a.second[i] - b.second[i]) > eps) return false;
    }
    return true;
  }


  bool compare_int_Coord_grid_pairs (const std::pair<int, Coord_grid > &a, const std::pair<int, Coord_grid > &b)
  {
    if (a.first < b.first) return true;
    if (a.first > b.first) return false;
    for (int i = 0; i < 3; i++) {
      if (a.second[i] > b.second[i]) return false;
      if (a.second[i] < b.second[i]) return true;
    }
    return false;
  }

  bool int_Coord_grid_pairs_equal (const std::pair<int, Coord_grid > &a, const std::pair<int, Coord_grid > &b)
  {

    if (a.first != b.first) return false;
    for (int i=0; i < 3; i++) {
      if (a.second[i] != b.second[i]) return false;
    }
    return true;
  }

  template <class T = Coord_frac> T cell_shift (const T& coord)
  {
    T ret;
    for (size_t i = 0; i < 3; i++) {
      ret[i] = floor(coord[i]);
    }
    return ret;
  }

  Coord_grid cell_shift (const Coord_grid& coord, const Grid_sampling& grid)
  {
    Coord_grid ret;
    ret = coord - coord.unit(grid);
    return ret;

  }
  }
%}

%ignore clipper::Symops::size() const;

// *INDENT-OFF*

%apply (double IN_ARRAY1[ANY]) {(double box_origin_xyz[3])};
%apply (double IN_ARRAY1[ANY]) {(double box_size_xyz[3])};
%apply (double IN_ARRAY1[ANY]) {(double box_res_xyz[3])};
%apply (long IN_ARRAY1[ANY]) {(long box_size_uvw[3])};
%apply (double IN_ARRAY2[ANY][ANY]) {(double model_bounds_xyz[8][3])};
%apply (double IN_ARRAY1[ANY]) {(double sides_xyz[3])};

// *INDENT-ON*

%feature("docstring") clipper::Symops "
  Stores a list of rotation/translation operators, which can be retrieved as
  Clipper RTOP_frac objects or numpy arrays of transformation matrices.
  "
%inline %{

  namespace clipper
  {
  /*! Convenience class to hold a vector of Symop objects, with
   *  Python getter/setter etc. functions.
   */

  class Symops
  {
    public:
      Symops()
      {
        size_ = 0;
      }
      Symops(std::vector<clipper::RTop_frac > ops)
      {
        symops_ = ops;
        size_ = symops_.size();
      }
      Symops(std::vector<clipper::Symop > ops)
      {
        size_ = 0;
        for (std::vector<clipper::Symop >::iterator it = ops.begin();
        it != ops.end(); ++it) {
          symops_.push_back(RTop_frac(it->rot(), it->trn()));
          unit_translations_.push_back(cell_shift<>(Coord_frac(it->trn())));
          size_ += 1;
        }
      }

      ~Symops() {}

      //! get element
      inline const clipper::RTop_frac operator []( const int& i ) const { return symops_.at(i); }

      const size_t& size() const {
        return size_;
      }

      clipper::RTop_frac __getitem__(int i)
      {
        i = (i < 0) ? size_ + i : i;
        if (i >= size_ || i < 0) {
          throw std::out_of_range("");
        }
        return symops_[i];
      }
      void __setitem__(int i, clipper::RTop_frac op)
      {
        i = (i < 0) ? size_ + i : i;
        if (i >= size_ || i < 0) {
          throw std::out_of_range("");
        }
        symops_[i] = op;
        unit_translations_[i] = cell_shift<>(Coord_frac(op.trn()));
      }
      size_t __len__()
      {
        return size_;
      }
      void append(clipper::RTop_frac op)
      {
        symops_.push_back(op);
        unit_translations_.push_back(cell_shift<>(Coord_frac(op.trn())));
        size_ += 1;
      }
      clipper::RTop_frac pop(int i)
      {
        clipper::RTop_frac ret = symops_[i];
        symops_.erase(symops_.begin() + i);
        unit_translations_.erase(unit_translations_.begin() + i);
        size_ -= 1;
        return ret;
      }
      Coord_frac cell_trans(int i)
      {
        return unit_translations_[i];
      }

      //! Return all fractional symops as a single array of 4x4 matrices
      void _all_matrices44_frac(double* numpy_array, int n1, int n2, int n3)
      {
        if (size_ != n1) {
          std::string errstring = "Target array length of " + std::to_string(n1) +
                   "does not match Symops array length of " + std::to_string(size_);
          throw std::length_error(errstring);
        }
        if (n2 !=4 || n3 != 4) {
          std::string errstring = "Target should be an nx4x4 numpy array!";
          throw std::length_error(errstring);
        }
        size_t count = 0;
        for (size_t i = 0; i < n1; i++) {
          RTop_frac thisop = symops_[i];
          for (size_t j = 0; j < n2; j++) {
            for (size_t k = 0; k < n3; k++, count++) {
              double thisval;
              // Fill in the bottom row
              if (j == 3) {
                if (k < 3) {
                  numpy_array[count] = 0;
                } else {
                  numpy_array[count] = 1;
                }
                continue;
              }
              // fill in the translations
              if (k == 3) {
                thisval = thisop.trn()[j];
              } else {
                // fill in the rotations
                thisval = thisop.rot()(j, k);
              }
              numpy_array[count] = thisval;
            }
          }
        }
      }

      //! Return all fractional symops as a single array of 3x4 matrices
      void _all_matrices34_frac(double* numpy_array, int n1, int n2, int n3)
      {
        if (size_ != n1) {
          std::string errstring = "Target array length of " + std::to_string(n1) +
                   "does not match Symops array length of " + std::to_string(size_);
          throw std::length_error(errstring);
        }
        if (n2 !=3 || n3 != 4) {
          std::string errstring = "Target should be an nx3x4 numpy array!";
          throw std::length_error(errstring);
        }
        size_t count = 0;
        for (size_t i = 0; i < n1; i++) {
          RTop_frac thisop = symops_[i];
          for (size_t j = 0; j < n2; j++) {
            for (size_t k = 0; k < n3; k++, count++) {
              double thisval;
              // fill in the translations
              if (k == 3) {
                thisval = thisop.trn()[j];
              } else {
                // fill in the rotations
                thisval = thisop.rot()(j, k);
              }
              numpy_array[count] = thisval;
            }
          }
        }
      }

      //! Return all orthographic symops as a single array of 4x4 matrices
      void _all_matrices44_orth(const Cell& cell, double* numpy_array, int n1, int n2, int n3)
      {
        if (size_ != n1) {
          std::string errstring = "Target array length of " + std::to_string(n1) +
                   "does not match Symops array length of " + std::to_string(size_);
          throw std::length_error(errstring);
        }
        if (n2 !=4 || n3 != 4) {
          std::string errstring = "Target should be an nx4x4 numpy array!";
          throw std::length_error(errstring);
        }
        size_t count = 0;
        for (size_t i = 0; i < n1; i++) {
          RTop_orth thisop = symops_[i].rtop_orth(cell);
          for (size_t j = 0; j < n2; j++) {
            for (size_t k = 0; k < n3; k++, count++) {
              double thisval;
              // Fill in the bottom row
              if (j == 3) {
                if (k < 3) {
                  numpy_array[count] = 0;
                } else {
                  numpy_array[count] = 1;
                }
                continue;
              }
              // fill in the translations
              if (k == 3) {
                thisval = thisop.trn()[j];
              } else {
                // fill in the rotations
                thisval = thisop.rot()(j, k);
              }
              numpy_array[count] = thisval;
            }
          }
        }
      }

      //! Return all orthographic symops as a single array of 3x4 matrices
      void _all_matrices34_orth(const Cell& cell, double* numpy_array, int n1, int n2, int n3)
      {
        if (size_ != n1) {
          std::string errstring = "Target array length of " + std::to_string(n1) +
                   "does not match Symops array length of " + std::to_string(size_);
          throw std::length_error(errstring);
        }
        if (n2 !=3 || n3 != 4) {
          std::string errstring = "Target should be an nx3x4 numpy array!";
          throw std::length_error(errstring);
        }
        size_t count = 0;
        for (size_t i = 0; i < n1; i++) {
          RTop_orth thisop = symops_[i].rtop_orth(cell);
          for (size_t j = 0; j < n2; j++) {
            for (size_t k = 0; k < n3; k++, count++) {
              double thisval;
              // fill in the translations
              if (k == 3) {
                thisval = thisop.trn()[j];
              } else {
                // fill in the rotations
                thisval = thisop.rot()(j, k);
              }
              numpy_array[count] = thisval;
            }
          }
        }
      }



    private:
      std::vector<clipper::RTop_frac >  symops_;
      std::vector<Coord_frac> unit_translations_;
      size_t size_;
  }; // class Symops
  } // namespace clipper
%} //inline

namespace clipper
{
  %extend Symops
  {
    %pythoncode %{
      def all_matrices_frac(self, format = '3x4'):
        '''
        Returns a Numpy array containing all the current symmetry operators in
        fractional coordinates.
        Args:
          format (string):
            One of '3x4' (default) or '4x4'. If '3x4', each operator is
            returned as:
            [[rot00, rot01, rot02, trn_u]
             [rot10, rot11, rot12, trn_v]
             [rot20, rot21, rot22, trn_w]]
            If '4x4', the [0,0,0,1] row necessary to form a complete
            affine transformation matrix is added.
        '''
        import numpy
        n = len(self)
        if format == '3x4':
          ret = numpy.empty((n,3,4), numpy.double)
          self._all_matrices34_frac(ret)
        elif format == '4x4':
          ret = numpy.empty((n,4,4), numpy.double)
          self._all_matrices44_frac(ret)
        else:
          raise TypeError("Format must be one of '3x4' or '4x4'!")
        return ret

      def all_matrices_orth(self, cell, format = '3x4'):
        '''
        Returns a Numpy array containing all the current symmetry operators in
        fractional coordinates.
        Args:
          cell:
            The clipper Cell object defining your crystal cell parameters
          format (string):
            One of '3x4' (default) or '4x4'. If '3x4', each operator is
            returned as:
            [[rot00, rot01, rot02, trn_x]
             [rot10, rot11, rot12, trn_y]
             [rot20, rot21, rot22, trn_z]]
            If '4x4', the [0,0,0,1] row necessary to form a complete
            affine transformation matrix is added.
        '''
        import numpy
        n = len(self)
        if format == '3x4':
          ret = numpy.empty((n,3,4), numpy.double)
          self._all_matrices34_orth(cell, ret)
        elif format == '4x4':
          ret = numpy.empty((n,4,4), numpy.double)
          self._all_matrices44_orth(cell, ret)
        else:
          raise TypeError("Format must be one of '3x4' or '4x4'!")
        return ret
    %}
  } // extend Symops
} // namespace clipper


%{
  namespace clipper {
  /*! Similar to Symops, but for integerised symops. Initialise from
   *  a Symops and a Grid_sampling object. Very useful internally, but
   *  I'm not sure if there's a need for it in Python. For now, let's not
   * bother wrapping it.
   */
  class Isymops
  {
    public:
      Isymops()
      {
        size_ = 0;
      }
      Isymops(const Symops& ops, const Grid_sampling& grid)
      {
        for (size_t i = 0; i < ops.size(); i++) {
          const RTop_frac& this_op = ops[i];
          /* There is no constructor for an Isymop from an RTop_frac,
           * so we need to convert to a Symop first. This somewhat
           * annoyingly discards the unit cell translations, so we'll
           * have to add them back to construct the final operator.
           */
          Symop this_symop = Symop(this_op);
          Isymop this_isym(this_symop, grid);
          Coord_grid this_trn = Coord_frac(this_op.trn()).coord_grid(grid);
          this_isym.trn() = this_trn;
          symops_.push_back(this_isym);
          unit_translations_.push_back(cell_shift(this_trn, grid));
        }
        size_ = symops_.size();
      }

      ~Isymops() {}

      //! get element
      inline const clipper::Isymop& operator []( const int& i ) const { return symops_.at(i); }

      const size_t& size() const {
        return size_;
      }

      clipper::Isymop __getitem__(int i)
      {
        i = (i < 0) ? size_ + i : i;
        if (i >= size_ || i < 0) {
          throw std::out_of_range("");
        }
        return symops_[i];
      }
      void __setitem__(int i, clipper::Isymop op, const Grid_sampling& grid)
      {
        i = (i < 0) ? size_ + i : i;
        if (i >= size_ || i < 0) {
          throw std::out_of_range("");
        }
        symops_[i] = op;
        unit_translations_[i] = cell_shift(Coord_grid(op.trn()), grid);
      }
      size_t __len__()
      {
        return size_;
      }
      void append(clipper::Isymop op, const Grid_sampling& grid)
      {
        symops_.push_back(op);
        unit_translations_.push_back(cell_shift(Coord_grid(op.trn()), grid));
        size_ += 1;
      }
      clipper::Isymop pop(int i)
      {
        clipper::Isymop ret = symops_[i];
        symops_.erase(symops_.begin() + i);
        unit_translations_.erase(unit_translations_.begin() + i);
        size_ -= 1;
        return ret;
      }
      Coord_grid cell_trans(int i)
      {
        return unit_translations_[i];
      }
    private:
      std::vector<clipper::Isymop >  symops_;
      std::vector<Coord_grid> unit_translations_;
      size_t size_;
  };
  } // namespace clipper
%} // inline
