/* mtz_types.cpp: ccp4 data types for the clipper libraries */
//C Copyright (C) 2000-2004 Kevin Cowtan and University of York


//L   This code is distributed under the terms and conditions of the
//L   CCP4 Program Suite Licence Agreement as a CCP4 Library.
//L   A copy of the CCP4 licence can be obtained by writing to the
//L   CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
#include "mtz_types.h"


namespace clipper {

MTZcrystal::MTZcrystal( const String& xname, const String& pname, const Cell& cell ) : Cell( cell )
{
  xname_ = xname;
  pname_ = pname;
}

const String& MTZcrystal::crystal_name() const
{
  return xname_;
}

const String& MTZcrystal::project_name() const
{
  return pname_;
}

MTZdataset::MTZdataset( const String& dname, const ftype& wavel )
{
  dname_ = dname;
  wavel_ = wavel;
}

const String& MTZdataset::dataset_name() const
{
  return dname_;
}

const ftype& MTZdataset::wavelength() const
{
  return wavel_;
}


} // namespace clipper
