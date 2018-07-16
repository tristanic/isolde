//  $Id: mmdb_math_fft.h $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2005-2013.
//
//    This library is free software: you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License version 3, modified in accordance with the provisions
//    of the license to address the requirements of UK law.
//
//    You should have received a copy of the modified GNU Lesser
//    General Public License along with this library. If not, copies
//    may be downloaded from http://www.ccp4.ac.uk/ccp4license.php
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//  =================================================================
//
//    12.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :   FFT <interface>
//       ~~~~~~~~~
//  **** Functions:  mmdb::math::FFT
//       ~~~~~~~~~   mmdb::math::RealFFT
//                   mmdb::math::TwoFFT
//                   mmdb::math::Convolve
//                   mmdb::math::mConvolve
//
//  (C) E.Krissinel  2005-2013
//
//  =================================================================
//

#ifndef  __FFT__
#define  __FFT__

#include "mmdb_mattype.h"

namespace mmdb  {

  namespace math  {

    extern void FFT      ( rvector data, int nn, bool Forward=true );

    extern void RealFFT  ( rvector data, int n,  bool Forward=true );

    extern void TwoFFT   ( rvector data1, rvector data2,
                           rvector fft1,  rvector fft2, int n );

    extern void Convolve ( rvector data, int n, rvector respns, int m,
                           rvector ans,  bool Conv=true );


    //   mConvolve ( data,n,m ) replaces array data[0..n-1] with the result
    // of m recursive convolutions (m>1) defined as
    //
    //      data_m = data (*) data_{m-1}
    //
    // where data_m is the result of mth convolution, data_0=data.
    // The definition of the convolution is
    //
    //      [a (*) b]_i = Sum_j { a_j * b_{i-j} }
    //
    //   On input, data[] is considered as containing the signal
    // sampled at both positive and negative times in the wrap-around
    // order, that is
    //
    //    data[i], 0<=i<n/2   signal sampled at times  dt*i
    //    data[i], n/2<=i<n   signal sampled at times -dt*(n-i)
    //
    // and the same wrap-around order is used to interprete the output
    // data. This means that if only m positive sampling times are
    // used, the length of data must be at least n=2*m, the rest being
    // padded with zeroes.
    //
    //   The number of sampling nodes n *must* be an integer power of
    // two, i.e. 2,4,8,16 ... .
    //
    extern void mConvolve ( rvector data, int n, int m );

  }  // namespace math

}  // namespace mmdb

#endif
