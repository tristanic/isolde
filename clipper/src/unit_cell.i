
%feature("docstring") clipper::Unit_Cell
  "Condensed description of a single unit cell based on an atomic model.\n
  Contains a reference fractional coordinate (typically corresponding
  to the centroid of the modelled asymmetric unit), a list of symmetry
  operators with the necessary translations added to fully pack a
  single unit cell, and a list of pre-computed inverse symmetry operators
  for convenience. Also provides methods to quickly find all the symmetry
  operators necessary to pack an arbitrary box in space with copies of
  your atomic model.
  "

%feature("docstring") clipper::Unit_Cell::symops
  "Returns a Symops object listing all the transformations mapping the
  reference model to other positions in the unit cell. "

%feature("docstring") clipper::Unit_Cell::inv_symops
  "Returns a Symops object listing all the transformations mapping each
  copy in the unit cell back to the reference model. "

%feature("docstring") clipper::Unit_Cell::ref_box
  "Returns a Grid_range object defining the smallest rhomboid in integer
  grid coordinates that fully encloses the atomic model. "

%feature("docstring") clipper::Unit_Cell::all_symops_in_box
  "all_symops_in_box(box_origin_xyz, box_size_uvw,
                    always_include_identity = False,
                    sample_frequency = 2)
      ---> Symops
  Returns a Symops object providing all symmetry operators which place
  some part of the atomic model within a box in grid space.
  NOTE: this routine is optimised for speed rather than strictness, and
  may occasionally return extra operators which place the model just
  outside the target box. Think of the result as a shortlist of operators
  upon which more detailed atom-by-atom searches can be done if necessary.
  Args:
    box_origin_xyz(numpy 1x3 double):
      Position of the bottom corner of the search box (in Angstroms)
    box_size_uvw(numpy 1x3 int):
      Box side lengths in grid coordinates
    always_include_identity(bool):
      If true, the identity operator will always be returned
    sample_frequency(int):
      The box will be searched in steps equal to the minimum side length
      of the reference box divided by this number. In almost all cases
      this should be left as is."

%rename ("_symops_ptr") clipper::Unit_Cell_clipper::Symops::symops();
%rename("_ref_asu_ptr") clipper::Unit_Cell::ref_asu();

%ignore clipper::Unit_Cell::isymops() const;
%ignore clipper::Unit_Cell::isymops();
%ignore clipper::Unit_Cell::inv_isymops() const;
%ignore clipper::Unit_Cell::inv_isymops();
%ignore clipper::Unit_Cell::ref_box_min_side() const;
%ignore clipper::Unit_Cell::find_symops_of_coord;
%ignore clipper::Unit_Cell::symops();
%ignore clipper::Unit_Cell::inv_symops();

%inline %{
  namespace clipper
  {
  //! Condensed unit cell description
  /*! Contains a reference fractional coordinate (e.g. corresponding
   *  to the centroid of the modelled asu), a list of symops
   *  with the necessary translations added to ensure all results
   *  cluster in the same unit cell, and a list of pre-computed inverse
   *  symops for convenience. Everything is defined relative to the
   *  given reference coordinate rather than Clipper's internal reference
   *  asu.
   */
  class Unit_Cell
  {
  public:
    Unit_Cell() {};
    ~Unit_Cell() {};

    //! Find the symops (with translations) necessary to pack a unit cell
    /*! The reference coordinate would typically be the centroid of
     *  the atomic model. Additionally, create a Grid_range object
     *  defining the grid that just encloses all atoms in the asymmetric
     *  unit.
     */
    Unit_Cell(Coord_frac ref, const Atom_list& atoms, const Cell& cell, const Spacegroup& sg, const Grid_sampling& grid, int padding=0)
    {
      ref_ = ref;
      grid_ = grid;
      cell_ = cell;
      sg_ = sg;
      /* Store the nearest grid coordinate to the input reference coordinate.
       * Make the equivalent grid, fractional and map reference coords
       * and store for repeated use.
       */
      grid_ref_ = ref.coord_grid(grid);
      cell_origin_ = grid_ref_ - grid_ref_.unit(grid);
      ref_ = grid_ref_.coord_frac(grid);

      /* Find the minimum and maximum grid coordinates of a box encompassing the atoms,
       * pad it by padding grid units either side, and make a Grid_range object from it.
       */
      Coord_grid ref_min = atoms[0].coord_orth().coord_frac(cell).coord_grid(grid);
      Coord_grid ref_max = ref_min;
      //Coord_grid pad(1,1,1);
      Coord_grid pad(padding,padding,padding);
      for (Atom_list::const_iterator it = atoms.begin(); it != atoms.end(); ++it) {
        const Coord_grid& thiscoord = it->coord_orth().coord_frac(cell).coord_grid(grid);
        for (size_t i = 0; i < 3; i++) {
          if (thiscoord[i] < ref_min[i]) ref_min[i] = thiscoord[i];
          else if (thiscoord[i] > ref_max[i]) ref_max[i] = thiscoord[i];
        }
      }
      //std::cerr << "Reference box corners: " << ref_min.format().c_str() << " " << ref_max.format().c_str() << std::endl;
      //std::cerr << ref_min.coord_frac(grid).coord_orth(cell).format().c_str() << ref_max.coord_frac(grid).coord_orth(cell).format().c_str() << std::endl;

      reference_model_bounds_ = Grid_range( ref_min-pad, ref_max+pad );
      //std::cerr << reference_model_bounds_.min().format().c_str() << reference_model_bounds_.max().format().c_str() << std::endl;

      // In order to catch all corner cases, we'll need to search the
      // block of 9 unit cells surrounding our reference model
      Coord_grid this_offset;
      for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
          for (int k = -1; k < 2; k++) {
            this_offset = Coord_grid(grid.nu()*i, grid.nv()*j, grid.nw()*k);
            ref_model_cell_origins_.push_back(cell_origin_ + this_offset);
          }
        }
      }
      /* Find the side of the reference box covering the shortest
       * distance. When searching a box for relevant symops, we'll
       * sample it in steps of half this size.
       */

      size_t nu, nv, nw;
      nu = reference_model_bounds_.nu();
      nv = reference_model_bounds_.nv();
      nw = reference_model_bounds_.nw();
      min_side_ = nv < nu ? (nw < nv ? nw : nv) : nu;

      for (size_t i = 0; i < sg.num_symops(); i++) {
        clipper::Symop this_symop = sg.symop(i);
        clipper::Coord_frac tc = ref_.transform(this_symop);
        // now we want to find the [int, int, int] shift that brings the
        // transformed coordinates back into the same unit cell as the reference
        // coordinate.
        clipper::RTop_frac new_op(this_symop.rot(), this_symop.trn());

        Coord_frac shift = tc - tc.lattice_copy_unit() - cell_origin_.coord_frac(grid);
        new_op.trn() -= shift;

        symops_.append(new_op);
      }
      for (size_t i = 0; i < symops_.size(); i++) {
        inv_symops_.append(symops_[i].inverse());
      }
      isymops_ = Isymops(symops_, grid);
      inv_isymops_ = Isymops(inv_symops_, grid);

    }
    // Symmetry operators mapping the reference coordinate to each
    // of the other asymmetric units in this unit cell
    const Symops& symops() const { return symops_;} //getter
    Symops& symops()
    {
      return symops_; //setter
    }
    // Symmetry operators mapping asymmetric units back to the
    // reference asu.
    const Symops& inv_symops() const {return inv_symops_;} //getter
    Symops& inv_symops()
    {
      return inv_symops_; //setter
    }
    // Reference coordinate (by default, the centroid of the atomic
    // model)
    const Isymops& isymops() const { return isymops_;} //getter
    Isymops& isymops()
    {
      return isymops_; //setter
    }
    const Isymops& inv_isymops() const {return inv_isymops_;} //getter
    Isymops& inv_isymops()
    {
      return inv_isymops_; //setter
    }

    const Coord_frac& ref() const {return ref_;} //getter
    Coord_frac& ref()
    {
      return ref_; //setter
    }
    // Asymmetric unit containing the reference coordinate
    const Grid_range& ref_box() const {return reference_model_bounds_;} //getter
    Grid_range& ref_box()
    {
      return reference_model_bounds_; //setter
    }

    // The length of the shortest side of the reference box
    const size_t& ref_box_min_side() const {return min_side_;}

    Coord_grid min() {return cell_origin_;}
    Coord_grid max() {return cell_origin_
            + Coord_grid(grid_.nu()-1, grid_.nv()-1, grid_.nw()-1);}



    //! Find all symops relevant to a given grid coordinate
    /*! Actually finds all inverse symops that map the coordinate back
     *  into our reference box, and appends a std::pair object containing
     *  the symop index and the unit cell translations to the input
     *  pairlist.
     *  Args:
     *    pairlist: vector to append the found symops to
     *    thecoord: input grid coordinate
     */
    void find_symops_of_coord(std::vector<std::pair<int, Coord_grid> >& pairlist, const Coord_grid& thecoord)
    {
      Coord_grid shift;
      Coord_grid t_coord;

      for (std::vector<Coord_grid>::iterator it = ref_model_cell_origins_.begin();
            it != ref_model_cell_origins_.end(); ++it) {

        // First work out how many unit cells we are away from our reference
        shift = thecoord - thecoord.unit(grid_) - *it;
        /* Now try the inverse symops, and work out which (if any) put us
         * within the bounds of our reference box. NOTE: since our box
         * is the limits of the model (not an ASU), it's entirely possible
         * for the coordinate to be in more than one box - or, indeed, no box
         * at all. Therefore, we should err on the side of caution and test
         * all symops for each point.
        */
        for (size_t i = 0; i < inv_symops_.size(); i++) {
          t_coord = (thecoord.unit(grid_) + *it).transform(inv_isymops_[i]);
          if (reference_model_bounds_.in_grid(t_coord))
            pairlist.push_back(std::pair<int, Coord_grid>(i, shift));
        }
      }
    }

    //! Find all symops mapping the reference asu to positions in a box
    /*! The box is defined by its origin in xyz coordinates, and the
     *  number of grid coordinates along each axis.
     */
    Symops all_symops_in_box(double box_origin_xyz[3], long box_size_uvw[3], bool always_include_identity = false, int sample_frequency = 2)
    {
      Symops ret;

      /* We'll do all the work of finding the symops in integer grid
       * coordinates for speed.
       */

      std::vector<std::pair<int, Coord_grid> > ops_and_translations;
      if (always_include_identity) {
        std::pair<int, Coord_grid> identity;
        identity.first = 0;
        identity.second = Coord_grid(0,0,0);
        ops_and_translations.push_back(identity);
      }

      Coord_grid ref_grid = ref_.coord_grid(grid_);
      int num_symops = isymops_.size();

      Coord_grid box_origin = Coord_orth(box_origin_xyz[0], box_origin_xyz[1], box_origin_xyz[2]).coord_frac(cell_).coord_grid(grid_);
      Coord_grid box_max = box_origin + Coord_grid(box_size_uvw[0], box_size_uvw[1], box_size_uvw[2]);

      /* We'll sample the box in steps equal to 1/2 the length of the
       * shortest side of the box encompassing the atomic model, making
       * sure we also capture the faces and corners.
       */

      size_t step_size = min_side_/sample_frequency;
      step_size = step_size == 0 ? 1 : step_size;
      bool u_done = false, v_done = false, w_done = false;
      int u,v,w;
      Coord_grid thiscg;
      u = box_origin[0];
      while (!u_done) {
        v=box_origin[1];
        v_done = false;
        if (u == box_max[0]) u_done = true;
        while (!v_done) {
          w=box_origin[2];
          w_done = false;
          if (v == box_max[1]) v_done = true;
          while (!w_done) {
            if (w == box_max[2]) w_done = true;
            thiscg = Coord_grid(u,v,w);
            this -> find_symops_of_coord(ops_and_translations, thiscg);
            w += step_size;
            if ( w > box_max[2] ) w = box_max[2];
          }
          v += step_size;
          if ( v > box_max[1] ) v = box_max[1];
        }
        u += step_size;
        if ( u > box_max[0] ) u = box_max[0];
      }
      // Sort the two vectors, and remove any duplicates
      std::sort(ops_and_translations.begin(), ops_and_translations.end(), compare_int_Coord_grid_pairs);

      ops_and_translations.erase( std::unique( ops_and_translations.begin(),
                                  ops_and_translations.end(),
                                  int_Coord_grid_pairs_equal
                                             ), ops_and_translations.end() );

      for (std::vector<std::pair<int, Coord_grid > >::iterator it = ops_and_translations.begin();
           it != ops_and_translations.end(); ++it) {
        RTop_frac thisop = symops_[it->first];
        thisop.trn() += (it->second).coord_frac(grid_);
        ret.append(thisop);

      }
      return ret;

    } // all_symops_in_box


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
    size_t min_side_; // The length of the smallest side of the grid box in reference_model_bounds
    /* Because life is never easy, the reference model will almost certainly
     * span across multiple unit cells (unlike the nice, neat map reference
     * asu used by Clipper. Therefore, to find symmetry equivalents we
     * need to test the main unit cell plus all 26 direct neighbours.
     */
    std::vector<Coord_grid> ref_model_cell_origins_;
  }; // Unit_Cell

  } // namespace clipper
%} // %inline


#ifdef PYTHON_PROPERTIES
namespace clipper
{
  %extend Unit_Cell
  {
    %pythoncode %{
      symops = property(symops)
      inv_symops = property(inv_symops)
      ref = property(ref)
      ref_box = property(ref_box)
      min = property(min)
      max = property(max)
    %}
  } // extend Symops
} // namespace clipper
#endif
