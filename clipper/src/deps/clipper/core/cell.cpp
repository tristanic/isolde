/* cell.cpp: class file for unit cell class */
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
/* This code is derived from the 'dm' source code */


#include "cell.h"


namespace clipper {


/*!
  \param a A axis in Angstroms.
  \param b B axis in Angstroms.
  \param c C axis in Angstroms.
  \param alpha Angle between B and C axes in radians or degrees, default=90
  \param beta  Angle between A and C axes in radians or degrees, default=90
  \param gamma Angle between A and C axes in radians or degrees, default=90 */
Cell_descr::Cell_descr( const ftype& a, const ftype& b, const ftype& c,
			const ftype& alpha, const ftype& beta,
			const ftype& gamma )
{
  a_=a; b_=b; c_=c; alpha_=alpha; beta_=beta; gamma_=gamma;
  if (alpha_>M_PI) alpha_=Util::d2rad(alpha_);
  if (beta_ >M_PI) beta_ =Util::d2rad(beta_ );
  if (gamma_>M_PI) gamma_=Util::d2rad(gamma_);
}

/*! \return The cell angle in degrees */
ftype Cell_descr::alpha_deg() const
{ return Util::rad2d(alpha_); }

/*! \return The cell angle in degrees */
ftype Cell_descr::beta_deg() const
{ return Util::rad2d(beta_ ); }

/*! \return The cell angle in degrees */
ftype Cell_descr::gamma_deg() const
{ return Util::rad2d(gamma_); }

/*! \return A string describing the cell */
String Cell_descr::format() const
{ return " Cell ("+String(a())+","+String(b())+","+String(c())+","+String(alpha_deg())+","+String(beta_deg())+","+String(gamma_deg())+")"; }


/*! Construct and initialise a metric tensor, given a set of real or
  reciprocal cell parameters.
  \param a Length of \b a axis in Angstroms or reciprocal Angstroms.
  \param b Length of \b b axis in Angstroms or reciprocal Angstroms.
  \param c Length of \b c axis in Angstroms or reciprocal Angstroms.
  \param alph Angle between \b b and \b c in radians.
  \param beta Angle between \b a and \b c in radians.
  \param gamm Angle between \b a and \b b in radians.
 */
Metric_tensor::Metric_tensor( const ftype& a, const ftype& b, const ftype& c, const ftype& alph, const ftype& beta, const ftype& gamm )
{
  m00 = a*a;
  m11 = b*b;
  m22 = c*c;
  m01 = 2.0*a*b*cos(gamm);
  m02 = 2.0*a*c*cos(beta);
  m12 = 2.0*b*c*cos(alph);
}

String Metric_tensor::format() const
{
  return "m00=" + String(m00) + " m11=" + String(m11) + " m22=" + String(m22) + " m01=" + String(m01) + " m02=" + String(m02) + " m12=" + String(m12);
}


// Cell object

/*! Initialise the Cell object from a cell description.
  \param cell_ The cell descirption. */
void Cell::init(const Cell_descr& cell_)
{
  // store cell (allows a/b/c/alpha/beta/gamma methods)
  a_ = cell_.a();
  b_ = cell_.b();
  c_ = cell_.c();
  alpha_ = cell_.alpha();
  beta_  = cell_.beta();
  gamma_ = cell_.gamma();

  // calculate cell volume (allows a_star../alpha_star.. methods)
  vol = a() * b() * c() * sqrt( 2.0*cos(alpha())*cos(beta())*cos(gamma())
				- cos(alpha())*cos(alpha())
				- cos( beta())*cos( beta())
				- cos(gamma())*cos(gamma()) + 1.0 );

  // deal with null cell
  if ( is_null() ) return;

  // othogonalisation+fractionisation matrices
  orthmat = Mat33<>::identity();
  orthmat(0,0) =  a();
  orthmat(0,1) =  b()*cos(gamma());
  orthmat(0,2) =  c()*cos( beta());
  orthmat(1,1) =  b()*sin(gamma());
  orthmat(1,2) = -c()*sin( beta())*cos(alpha_star());
  orthmat(2,2) =  c()*sin( beta())*sin(alpha_star());
  //std::cout << sin( beta()) << " " << sin(alpha_star()) << " " << alpha_star() << "\n";
  fracmat = orthmat.inverse();

  // metric tensors
  realmetric = Metric_tensor( a(), b(), c(),
			      alpha(), beta(), gamma() );
  recimetric = Metric_tensor( a_star(), b_star(), c_star(),
			      alpha_star(), beta_star(), gamma_star() );
}

/*! \return true if the object has not been initalised. */
bool Cell::is_null() const
{ return ( vol <= 0.0 ); }

ftype Cell::a_star() const
{ return b()*c()*sin(alpha())/vol; }

ftype Cell::b_star() const
{ return c()*a()*sin(beta())/vol; }

ftype Cell::c_star() const
{ return a()*b()*sin(gamma())/vol; }

ftype Cell::alpha_star() const
{ return acos( (cos(gamma())*cos(beta())-cos(alpha())) /
	       (sin(beta())*sin(gamma())) ); }

ftype Cell::beta_star() const
{ return acos( (cos(alpha())*cos(gamma())-cos(beta())) /
	       (sin(gamma())*sin(alpha())) ); }

ftype Cell::gamma_star() const
{ return acos( (cos(beta())*cos(alpha())-cos(gamma())) /
	       (sin(alpha())*sin(beta())) ); }

/*! Two cells disagree if the difference in their orthogonalisation
  matrices is sufficient to map a reflection from one cell onto a
  different reflection in the other cell at the given tolerance, which
  is the resolution of the reflection in Angstroms.
  \param other The other cell to compare.
  \param tol The tolerance, in Angstroms. */
bool Cell::equals( const Cell& other, const ftype tol ) const
{
  if ( is_null() || other.is_null() ) return false;
  ftype s = 0.0;
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      s += pow( fracmat(i,j) - other.fracmat(i,j), 2.0 );
  return s < ( pow( tol, 2.0 ) / pow( vol, 1.333 ) );
}

void Cell::debug() const
{
  std::cout << descr().format() << "\n";
  std::cout << "Vol" << vol << "\n";
  std::cout << "Orth mat\n" << orthmat.format() << "\n";
  std::cout << "Frac mat\n" << fracmat.format() << "\n";
  std::cout << "Prod mat\n" << (orthmat*fracmat).format() << "\n";
  std::cout << "Real metric " << realmetric.format() << "\n";
  std::cout << "Reci metric " << recimetric.format() << "\n";
  //std::cout << "Real metric\n"+(orthmat.transpose()*orthmat).format()+"\n";
  //std::cout << "Reci metric\n"+(fracmat*fracmat.transpose()).format()+"\n";
}


} // namespace clipper
