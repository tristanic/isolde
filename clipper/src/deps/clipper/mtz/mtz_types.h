/*! \file lib/mtz_types.h
    Header file for CCP4 data types for the clipper libraries
*/
//C Copyright (C) 2000-2004 Kevin Cowtan and University of York


//L   This code is distributed under the terms and conditions of the
//L   CCP4 Program Suite Licence Agreement as a CCP4 Library.
//L   A copy of the CCP4 licence can be obtained by writing to the
//L   CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
#ifndef CLIPPER_MTZ_TYPES
#define CLIPPER_MTZ_TYPES


#include "../core/container_hkl.h"


namespace clipper
{

  //! CCP4 crystal object
  class MTZcrystal : public Cell
  {
  public:
    //! null constructor
    MTZcrystal() {}
    //! constructor: takes crystal+project name and cell
    MTZcrystal( const String& xname, const String& pname, const Cell& cell );
    //! get crystal name
    const String& crystal_name() const;
    //! get project name
    const String& project_name() const;
  protected:
    String xname_, pname_;
  };

  //! CCP4 dataset object
  class MTZdataset
  {
  public:
    //! null constructor
    MTZdataset() {}
    //! constructor: takes wavelength
    MTZdataset( const String& dname, const ftype& wavel );
    //! get dataset name
    const String& dataset_name() const;
    //! get wavelength
    const ftype& wavelength() const;
  protected:
    String dname_;
    ftype wavel_;
  };

  //! CMTZcrystal identifier
  /*!
    CMTZcrystal: This has a name and a cell.
    It overrides the base cell for any HKL_datas below it, and mirrors the
    mtz++ crystal element.
  */
  class CMTZcrystal : public Container, public MTZcrystal
  {
  public:
    //! null constructor
    CMTZcrystal() {}
    //! constructor: from MTZcrystal
    CMTZcrystal( Container& parent, const String& name, const MTZcrystal& xtl ) : Container( parent, name ), MTZcrystal( xtl ) {}
  };

  //! CMTZdataset identifier
  /*!
    CMTZdataset: This has a name and a wavelength.
    It gives the wavelength for any HKL_datas below it, and mirrors the
    mtz++ dataset element.
  */
  class CMTZdataset : public Container, public MTZdataset
  {
  public: 
    //! null constructor
    CMTZdataset() {}
    //! constructor: from MTZdataset
    CMTZdataset( Container& parent, const String& name, const MTZdataset& set ) : Container( parent, name ), MTZdataset( set ) {}
  };

} // namespace clipper

#endif
