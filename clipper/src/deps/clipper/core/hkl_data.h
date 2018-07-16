/*! \file lib/hkl_data.h
    Header file for reflection data class
*/
//C Copyright (C) 2000-2006 Kevin Cowtan and University of York
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA


#ifndef CLIPPER_HKL_DATA
#define CLIPPER_HKL_DATA


#include "hkl_info.h"
#include "../imex.h"

namespace clipper
{

  class CLIPPER_IMEX HKL_data_cacheobj : public HKL_info
  {
  public:
    class Key
    {
    public:
      Key( const Spgr_descr& spgr_descr, const Cell& cell_descr, const HKL_sampling& hkl_sam ) : spgr_descr_(spgr_descr), cell_descr_(cell_descr), hkl_sampling_(hkl_sam) {}
      const Spgr_descr& spgr_descr() const { return spgr_descr_; }
      const Cell_descr& cell_descr() const { return cell_descr_; }
      const HKL_sampling& hkl_sampling() const { return hkl_sampling_; }
    private:
      Spgr_descr spgr_descr_;
      Cell_descr cell_descr_;
      HKL_sampling hkl_sampling_;
    };

    HKL_data_cacheobj( const Key& hkl_data_cachekey );
    bool matches( const Key& hkl_data_cachekey ) const;
    String format() const;
    static Mutex mutex;                 //!< thread safety
  private:
    Key key;
  };


  //! Reflection data type objects.
  /*!
    A class from which data type objects are usually derived

    To define a new type for use in an HKL_data structure, subclass
    this class and override the following methods:
    - constructor   - initialises to NAN
    - type()        - returns type name, which is a list of the contained data
    - friedel()     - applies Fridel transformation
    - shift_phase() - applies phase shift transformation
    - missing()     - checks if data is present

    - data_size()   - number of data elements in this type
    - data_names()  - names of data elements in this type
    - data_export() - conversion to array (for I/O)
    - data_import() - conversion from array (for I/O)

    - scale()       - (OPTIONAL) apply magnitude scale factor to data

    \note polymorphism is NOT used here because virtual tables would be
    to expensive for every individual reflection, both in terms of
    space and cpu cycles.
  */
  class CLIPPER_IMEX Datatype_base
  {
  protected:
    //! initialise to 'missing'
    //-- Datatype_base();
    //! initialise to 'missing' (all elements are set to null)
    void set_null();
    //! return the name of this data type
    static String type();
    //! apply Friedel transformation to the data
    void friedel();
    //! apply phase shift to the data
    void shift_phase(const ftype& dphi);
    //! return true if data is missing
    bool missing() const;
    //! return number of data elements in this type
    static int data_size();
    //! return names of data elements in this type
    static String data_names();
    //! conversion to array (for I/O)
    void data_export( xtype array[] ) const;
    //! conversion from array (for I/O)
    void data_import( const xtype array[] );
  };


  //! HKL_data_base
  /*!
    This is the virtual base for the typed hkl_data objects. It exists
    to guarantee and interface by which data can be managed without
    knowledge of the specific data type.
  */

  class CLIPPER_IMEX HKL_data_base
  {
   public:
    // Coordinate reference types
    //! Basic HKL_reference_index: see HKL_info
    typedef HKL_info::HKL_reference_index HKL_reference_index;
    //! HKL HKL_reference_index: see HKL_info
    typedef HKL_info::HKL_reference_coord HKL_reference_coord;

    //! initialiser: from parent hkl_info and cell
    virtual void init( const HKL_info& hkl_info, const Cell& cell );
    //! initialiser: from another hkl_data
    virtual void init( const HKL_data_base& hkl_data );
    //! [CLIPPER2] initialiser: from spacegroup, cell, and HKL_sampling
    virtual void init( const Spacegroup& spacegroup, const Cell& cell, const HKL_sampling& hkl_sampling );

    // generic methods
    //! test if object has been initialised
    bool is_null() const;

    //! get the parent HKL_info object
    const HKL_info& base_hkl_info() const { return *parent_hkl_info; }
    //! get the parent cell
    const Cell& base_cell() const { return *parent_cell; }

    //! [CLIPPER2] get spacegroup
    const Spacegroup& spacegroup() const { return spacegroup_; }
    //! [CLIPPER2] get cell
    const Cell& cell() const { return cell_; }
    //! [CLIPPER2] get resolution
    const Resolution& resolution() const { return resolution_; }
    //! [CLIPPER2] get HKL_sampling
    const HKL_sampling& hkl_sampling() const { return hkl_sampling_; }
    //! [CLIPPER2] get HKL_info object
    const HKL_info& hkl_info() const { return *parent_hkl_info; }

    //! get resolution by reflection index (based on true cell)
    ftype invresolsq( const int& index ) const;
    //! get resolution limits of the list (based on true cell and missing data)
    Range<ftype> invresolsq_range() const;
    //! get number of observations in this list (based on missing data)
    int num_obs() const;

    //! update: synchornize info with parent HKL_info
    virtual void update() = 0;
    //! get data type (a list of names corresponding to the im/export values)
    virtual String type() const = 0;
    //! check if a data entry in the list is marked as 'missing'
    virtual bool missing(const int& index) const = 0;
    //! set data entry in the list to its null value
    virtual void set_null(const int& index) = 0;
    //! return number of data elements in this type
    virtual int data_size() const = 0;
    //! return names of data elements in this type
    virtual String data_names() const = 0;
    //! conversion to array (for I/O)
    virtual void data_export( const HKL& hkl, xtype array[] ) const = 0;
    //! conversion from array (for I/O)
    virtual void data_import( const HKL& hkl, const xtype array[] ) = 0;
    //! mask the data by marking any data missing in 'mask' as missing
    virtual void mask(const HKL_data_base& mask) = 0;

    //! return HKL_reference_index pointing to first reflection
    HKL_reference_index first() const;
    //! return HKL_reference_index pointing to first non-missing data
    HKL_reference_index first_data() const;
    //! increment HKL_reference_index to next non-missing data
    HKL_reference_index& next_data( HKL_reference_index& ih ) const;

    void debug() const;

   protected:
    const HKL_info* parent_hkl_info;
    const Cell* parent_cell;
    bool cell_matches_parent;

    // clipper2 members
    ObjectCache<HKL_data_cacheobj>::Reference cacheref;  //!< object cache ref
    Spacegroup spacegroup_;
    Cell cell_;
    HKL_sampling hkl_sampling_;
    Resolution resolution_;

    //! null constructor
    HKL_data_base();
    //! destructor
    virtual ~HKL_data_base() {}
  };


  //! HKL_data<>
  /*!
    An actual hkl_data object, containing actual data of type T. This
    implements the generic interface, and in addition provides
    type-specific access functions.

    \note The following methods are inherited from HKL_data_base but are documented here for convenience: base_hkl_info(), base_cell(), invresolsq(), invresolsq_range(), num_obs(), first(), first_data(), next_data().
  */
  template<class T> class HKL_data : public HKL_data_base
  {
  public:
    //! null constructor
    HKL_data() {}
    //! constructor: from parent hkl_info
    explicit HKL_data( const HKL_info& hkl_info );
    //! constructor: from parent hkl_info and cell
    HKL_data( const HKL_info& hkl_info, const Cell& cell );
    //! [CLIPPER2] constructor: from spacegroup, cell and hkl_sampling
    HKL_data( const Spacegroup& spacegroup, const Cell& cell, const HKL_sampling& hkl_sampling );
    //! [CLIPPER2] constructor: from another HKL_data object
    explicit HKL_data( const HKL_data_base& hkl_data );

    //! initialiser: from parent hkl_info and cell
    void init( const HKL_info& hkl_info, const Cell& cell );
    //! [CLIPPER2] initialiser: from spacegroup, cell, and HKL_sampling
    void init( const Spacegroup& spacegroup, const Cell& cell, const HKL_sampling& hkl_sampling );
    //! [CLIPPER2] initialiser: from another HKL_data object
    void init( const HKL_data_base& hkl_data );
    //! update: synchornize info with parent HKL_info
    void update();

    // type specific methods
    String type() const { return T::type(); }
    bool missing(const int& index) const { return list[index].missing(); }
    void set_null(const int& index) { list[index].set_null(); }
    int data_size() const { return T::data_size(); }
    String data_names() const { return T::data_names(); }
    void data_export( const HKL& hkl, xtype array[] ) const
      { T datum; get_data( hkl, datum ); datum.data_export( array ); }
    void data_import( const HKL& hkl, const xtype array[] )
      { T datum; datum.data_import( array ); set_data( hkl, datum ); }
    void mask(const HKL_data_base& mask);

    // data access methods: by HKL_reference_index
    //! get data by reflection HKL_reference_index
    const T& operator[] (const HKL_info::HKL_reference_index& i) const
      { return list[i.index()]; }
    //! set data by reflection HKL_reference_index
    T& operator[] (const HKL_info::HKL_reference_index& i)
      { return list[i.index()]; }

    // data access methods: by HKL_reference_coord
    //! get data by HKL_reference_coord
    T operator[] (const HKL_info::HKL_reference_coord& ih) const;
    //! get data by HKL_reference_coord (returns false if no equivalent hkl)
    bool get_data(const HKL_info::HKL_reference_coord& ih, T& data) const;
    //! set data by HKL_reference_coord (returns false if no equivalent hkl)
    bool set_data(const HKL_info::HKL_reference_coord& ih, const T& data);

    // data access methods: by index
    //! get data by reflection index
    const T& operator[] (const int& index) const { return list[index]; }
    //! set data by reflection index
    T& operator[] (const int& index) { return list[index]; }

    // data access methods: by hkl
    //! get data by hkl (returns missing if no equivalent hkl)
    T operator[] (const HKL& hkl) const;
    //! get data by hkl (returns false if no equivalent hkl)
    bool get_data(const HKL& hkl, T& data) const;
    //! set data by hkl (returns false if no equivalent hkl)
    bool set_data(const HKL& hkl, const T& data);

    // COMPUTATION OPERATORS
    //! Basic computation: fill this data list by function call
    template<class C> void compute( const C& op )
      { for (HKL_info::HKL_reference_index ih=parent_hkl_info->first(); !ih.last(); ih.next()) list[ih.index()] = op( ih ); }
    //! Unary computation: fill this data list by computation from another
    template<class S, class C> void compute( const HKL_data<S>& src, const C& op )
      { for (HKL_info::HKL_reference_index ih=parent_hkl_info->first(); !ih.last(); ih.next()) list[ih.index()] = op( ih, src[ih] ); }
    //! Binary computation: fill this data list by computation from another
    template<class S1, class S2, class C> void compute( const HKL_data<S1>& src1, const HKL_data<S2>& src2, const C& op )
      { for (HKL_info::HKL_reference_index ih=parent_hkl_info->first(); !ih.last(); ih.next()) list[ih.index()] = op( ih, src1[ih], src2[ih] ); }

    // inherited functions lists for documentation purposes
    //-- const HKL_info& base_hkl_info() const;
    //-- const Cell& base_cell() const;
    //-- const ftype invresolsq(const int& index) const;
    //-- const Range<ftype> invresolsq_range() const;
    //-- const int num_obs() const;
    //-- HKL_reference_index first() const;
    //-- HKL_reference_index first_data() const;
    //-- HKL_reference_index& next_data( HKL_reference_index& ih ) const;

    //! assignment operator: copies the data from another list
    HKL_data<T>& operator =( const HKL_data<T>& other );
    //! assignment operator: assigns a single value to the whole list
    HKL_data<T>& operator =( const T& value );

    void debug() const;

  protected:
    // members
    std::vector<T> list;
  };







  // Template implementations

  ftype HKL_info::HKL_reference_base::invresolsq( const HKL_data_base& hkldata ) const
    { return hkldata.invresolsq( index_ ); }

  /*! Construct the object using a given reflection list and cell.
    \param hkl_info The reflection list object. */
  template<class T> HKL_data<T>::HKL_data( const HKL_info& hkl_info )
  {
    init( hkl_info, hkl_info.cell() );
  }

  /*! Construct the object using a given reflection list and cell.
    \param hkl_info The reflection list object.
    \param cell The unit cell for this datalist. */
  template<class T> HKL_data<T>::HKL_data( const HKL_info& hkl_info, const Cell& cell )
  {
    init( hkl_info, cell );
  }

  /*! Construct the object using a given spacegroup, cell, and sampling.
  \param spacegroup The spacegroup for this datalist.
  \param cell The unit cell for this datalist.
  \param hkl_sampling The reflection list description. */
  template<class T> HKL_data<T>::HKL_data( const Spacegroup& spacegroup, const Cell& cell, const HKL_sampling& hkl_sampling )
  {
    init( spacegroup, cell, hkl_sampling );
  }

  /*! Construct the object using a given HKL_data object. The
    properties of the object (spacegroup, cell, sampling) are the
    copied, but the actual data is not.
  \param hkl_data The HKL_data object to provide the data. */
  template<class T> HKL_data<T>::HKL_data( const HKL_data_base& hkl_data )
  {
    init( hkl_data );
  }

  /*! Initialise the object using a given reflection list and cell.
    \param hkl_info The reflection list object.
    \param cell The unit cell for this datalist. */
  template<class T> void HKL_data<T>::init( const HKL_info& hkl_info, const Cell& cell )
  {
    HKL_data_base::init( hkl_info, cell );
    update();
  }

  /*! Initialise the object using a given spacegroup, cell, and sampling.
    \param spacegroup The spacegroup for this datalist.
    \param cell The unit cell for this datalist.
    \param hkl_sampling The reflection list description. */
  template<class T> void HKL_data<T>::init( const Spacegroup& spacegroup, const Cell& cell, const HKL_sampling& hkl_sampling )
  {
    HKL_data_base::init( spacegroup, cell, hkl_sampling );
    update();
  }

  /*! Initialise the object using a given HKL_data object. The
    properties of the object (spacegroup, cell, sampling) are the
    copied, but the actual data is not.
  \param hkl_data The HKL_data object to provide the data. */
  template<class T> void HKL_data<T>::init( const HKL_data_base& hkl_data )
  {
    HKL_data_base::init( hkl_data );
    update();
  }

  /*! The datalist is resized if necessary to match the parent. */
  template<class T> void HKL_data<T>::update()
  {
    if ( parent_hkl_info != NULL ) {
      T null; null.set_null();
      list.resize( parent_hkl_info->num_reflections(), null );
    }
  }

  /*! For each data element, if the corresponding element in \c mask
    is missing, then that element in this list is also set to
    missing.
    \param mask The list to provide the mask. */
  template<class T> void HKL_data<T>::mask(const HKL_data_base& mask)
  {
    T null; null.set_null();
    for ( unsigned int i = 0; i < list.size(); i++ )
      if ( mask.missing(i) ) list[i] = null;
  }

  /*! If a symmetry mate of the requested HKL exists in the list, then
    the correct symmetry transformations are applied and the data is
    returned, otherwise the value of 'missing' for the datatype is returned.
    \param ih The reference to the HKL.
    \return The data, or 'missing'. */
  template<class T> T HKL_data<T>::operator[] (const HKL_info::HKL_reference_coord& ih) const
  {
    if ( ih.index() < 0 ) { T null; null.set_null(); return null; }
    T data = list[ih.index()];
    if ( ih.friedel() ) data.friedel();
    data.shift_phase( -ih.hkl().sym_phase_shift( parent_hkl_info->spacegroup().symop(ih.sym()) ) );
    return data;
  }

  /*! If a symmetry mate of the requested HKL exists in the list, then
    the correct symmetry transformations are applied and the data is
    returned, otherwise the value of 'missing' for the datatype is returned.
    \param ih The reference to the HKL.
    \param data Returned with the value of the data.
    \return true if the data was returned. */
  template<class T> bool HKL_data<T>::get_data(const HKL_info::HKL_reference_coord& ih, T& data) const
  {
    if ( ih.index() < 0 ) { data.set_null(); return false; }
    data = list[ih.index()];
    if ( ih.friedel() ) data.friedel();
    data.shift_phase( -ih.hkl().sym_phase_shift( parent_hkl_info->spacegroup().symop(ih.sym()) ) );
    return true;
  }

  /*! If a symmetry mate of the requested HKL exists in the list, then
    the correct symmetry transformations are applied and data is set
    to the supplied values, otherwise the function returns false.
    \param ih The reference to the HKL.
    \param data Value of the data to set.
    \return true if the data was set. */
  template<class T> bool HKL_data<T>::set_data(const HKL_info::HKL_reference_coord& ih, const T& data)
  {
    if ( ih.index() < 0 ) return false;
    T& ldata = list[ih.index()];
    ldata = data;
    ldata.shift_phase( ih.hkl().sym_phase_shift( parent_hkl_info->spacegroup().symop(ih.sym()) ) );
    if ( ih.friedel() ) ldata.friedel();
    return true;
  }

  /*! If a symmetry mate of the requested HKL exists in the list, then
    the correct symmetry transformations are applied and the data is
    returned, otherwise the value of 'missing' for the datatype is returned.
    \param hkl The reflection HKL.
    \return The data, or 'missing'. */
  template<class T> T HKL_data<T>::operator[] (const HKL& hkl) const
  {
    int index, sym; bool friedel;

    index = parent_hkl_info->index_of( parent_hkl_info->
				       find_sym(hkl, sym, friedel) );
    if ( index < 0 ) { T null; null.set_null(); return null; }
    T data = list[index];
    if (friedel) data.friedel();
    data.shift_phase( -hkl.sym_phase_shift( parent_hkl_info->spacegroup().symop(sym) ) );
    return data;
  }

  /*! If a symmetry mate of the requested HKL exists in the list, then
    the correct symmetry transformations are applied and the supplied
    datatype is set, otherwise the function returns false.
    \param hkl The reflection HKL.
    \param data Returned with the value of the data.
    \return true if the data was returned. */
  template<class T> bool HKL_data<T>::get_data(const HKL& hkl, T& data) const
  {
    int index, sym; bool friedel;

    index = parent_hkl_info->index_of( parent_hkl_info->
					find_sym(hkl, sym, friedel) );
    if ( index < 0 ) { data.set_null(); return false; }
    data = list[index];
    if (friedel) data.friedel();
    data.shift_phase( -hkl.sym_phase_shift(parent_hkl_info->spacegroup().symop(sym)) );
    return true;
  }

  /*! If a symmetry mate of the requested HKL exists in the list, then
    the correct symmetry transformations are applied and data is set
    to the supplied values, otherwise the function returns false.
    \param hkl The reflection HKL.
    \param data Value of the data to set.
    \return true if the data was set. */
  template<class T> bool HKL_data<T>::set_data(const HKL& hkl, const T& data_)
  {
    int index, sym; bool friedel;
    index = parent_hkl_info->index_of( parent_hkl_info->
				       find_sym(hkl, sym, friedel) );
    if ( index < 0 ) { return false; }
    T& ldata = list[index];
    ldata = data_;
    ldata.shift_phase( hkl.sym_phase_shift(parent_hkl_info->spacegroup().symop(sym)) );
    if (friedel) ldata.friedel();
    return true;
  }


  /*! The data list is copied from the assignment source to the
    target. If the target does not have a defined HKL_info, then that
    and the Cell are copied as well. If however the target does have a
    defined HKL_info the HKL_info objects are compared, and if they do
    not match an exception is thrown.
    \param other The datalist to copy.
    \return This list. */
  template<class T> HKL_data<T>& HKL_data<T>::operator =( const HKL_data<T>& other )
  {
    if ( parent_hkl_info == NULL ) {
      init( other );
    } else {
      if ( parent_hkl_info != other.parent_hkl_info )
	Message::message( Message_fatal( "HKL_data<T>: mismatched parent HKL_info is assignment" ) );
    }
    list = other.list;
    return *this;
  }


  /*! All values, including missing values, are overwritten by the value.
    \param value The value to which the list is to be set. */
  template<class T> HKL_data<T>& HKL_data<T>::operator =( const T& value )
  {
    for ( unsigned int i = 0; i < list.size(); i++ ) list[i] = value;
    return *this;
  }


  template<class T> void HKL_data<T>::debug() const
  {
    HKL_data_base::debug();
    std::cout << "Size " << list.size() << "\n";
  }


} // namespace clipper

#endif
