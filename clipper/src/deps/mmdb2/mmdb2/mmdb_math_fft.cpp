//  $Id: mmdb_math_fft.cpp $
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
//  **** Module  :   FFT <implementation>
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

#include <math.h>

#include "mmdb_math_fft.h"

namespace mmdb  {

  namespace math  {

    void  FFT ( rvector data, int nn, bool Forward )  {
    //  Replaces data[1:2*nn] by its discrete Fourier transform,
    //  if Forward is true; or replaces data[1:2*nn] by nn times its
    //  inverse discrete Fourier transform if Forward is false.
    //    On call,
    //  data[i], i=1,3,5 ... nn-1 are real parts of function,
    //  data[i], i=2,4,6 ... nn   are imaginary parts.
    //    nn MUST be an integer power of 2 (this is not checked for!).
    //    User should allocate data with GetVectorMemory and deallocate
    //  it with FreeVectorMemory to assure correct managing of indices
    //  1..2*nn (not 0..2*nn-1).
    //    On return,
    //  data[i], i=1,3,5 ... nn-1 are real parts, and
    //  data[i], i=2,4,6 ... nn   are imaginary parts of positive part
    //           of spectra, with frequences
    //            0,1/(nn*D),2/(nn*D) ... (nn/2-1)/(nn*D);
    //  data[nn+1] and data[nn+2] are real and imaginary parts for
    //           the frequances +/-1/(2*D)
    //  data[i], i=nn+3,nn+5,nn+7 ... 2*nn-1  are real parts, and
    //  data[i], i=nn+4,nn+6,nn+8 ... 2*nn    are imaginary parts of
    //           negative part of the spectra with frequences
    //            -(nn/2-1)/(nn*D), -(nn/2-2)/(nn*D), -1/(nn*D)
    //
    int         i,istep,j,m,mmax,n;
    realtype    tempi,tempr;
    long double theta,wi,wpi,wpr,wr,wtemp;  // this should be of
                                            // maximal precision
      n = 2*nn;
      j = 1;
      for (i=1;i<=n;i+=2)  {
        if (j>i)  {
          tempr     = data[j];
          tempi     = data[j+1];
          data[j]   = data[i];
          data[j+1] = data[i+1];
          data[i]   = tempr;
          data[i+1] = tempi;
        }
        m = n/2;
        while ((m>=2) && (j>m)) {
          j -= m;
          m /= 2;
        }
        j += m;
      }
      mmax = 2;
      while (n>mmax)  {
        istep = 2*mmax;
        theta = 2.0*Pi/mmax;
        if (!Forward)  theta = -theta;
        wpr = sin(0.5*theta);
        wpr = -2.0*wpr*wpr;
        wpi = sin(theta);
        wr  = 1.0;
        wi  = 0.0;
        for (m=1;m<=mmax;m+=2)  {
          for (i=m;i<=n;i+=istep)  {
            j     = i + mmax;
            tempr = wr*data[j]   - wi*data[j+1];
            tempi = wr*data[j+1] + wi*data[j];
            data[j]   = data[i]   - tempr;
            data[j+1] = data[i+1] - tempi;
            data[i]   = data[i]   + tempr;
            data[i+1] = data[i+1] + tempi;
          }
          wtemp = wr;
          wr    = wr*wpr - wi*wpi + wr;
          wi    = wi*wpr + wtemp*wpi + wi;
        }
        mmax = istep;
      }
    }

    void  RealFFT ( rvector data, int n, bool Forward )  {
    //    Calculates the Fourier transform of a set of n real-valued data
    //  points. Replaces this data (which is stored in array data[1:n])
    //  by the positive frequency half of its complex Fourier transform.
    //  The real-valued first and last components of the complex transform
    //  are returned as elements data[1] and data[2], respectively.
    //  n MUST be a power of 2. This routine also calculates the
    //  inverse transform of a complex data array if it is the transform
    //  of real data (Result in this case must be multiplied by 2/n).
    //    Array data should be allocated with GetVectorMemory.
    //
    int         i,i1,i2,i3,i4,n2p3;
    realtype    c1,c2,h1i,h1r,h2i,h2r;
    long double theta,wi,wpi,wpr,wr,wtemp;

      theta = 2.0*Pi/n;
      c1    = 0.5;
      if (Forward) {
        c2 = -0.5;
        FFT ( data,n/2,true );
      } else  {
        c2    = 0.5;
        theta = -theta;
      }
      wpr  = sin(0.5*theta);
      wpr  = -2.0*wpr*wpr;
      wpi  = sin(theta);
      wr   = 1.0 + wpr;
      wi   = wpi;
      n2p3 = n + 3;
      for (i=2;i<=n/4;i++)  {
        i1 = 2*i - 1;
        i2 = i1  + 1;
        i3 = n2p3 - i2;
        i4 = i3 + 1;
        h1r = c1*(data[i1] + data[i3]);
        h1i = c1*(data[i2] - data[i4]);
        h2r = -c2*(data[i2] + data[i4]);
        h2i = c2*(data[i1] - data[i3]);
        data[i1] = h1r + wr*h2r - wi*h2i;
        data[i2] = h1i + wr*h2i + wi*h2r;
        data[i3] = h1r - wr*h2r + wi*h2i;
        data[i4] = -h1i + wr*h2i + wi*h2r;
        wtemp = wr;
        wr = wr*wpr - wi*wpi + wr;
        wi = wi*wpr + wtemp*wpi + wi;
      }
      if (Forward)  {
        h1r = data[1];
        data[1] = h1r + data[2];
        data[2] = h1r - data[2];
      } else  {
        h1r = data[1];
        data[1] = c1*(h1r+data[2]);
        data[2] = c1*(h1r-data[2]);
        FFT ( data,n/2,false );
      }
    }

    void TwoFFT ( rvector data1, rvector data2,
                  rvector fft1,  rvector fft2, int n )  {
    //  Given two real input arrays data1[1:n] and data2[1:n],
    // this routine calls FFT and returns two "complex" output
    // arrays fft1[1:2*n] and fft2[1:2*n] (2*i-1 ith real,
    // 2*i ith imaginary), which contain the discrete Fourier
    // transforms of the respective data arrays.  n MUST be
    // an integer power of 2.
    int      i,j,n2, bj,bn;
    realtype h1r,h1i,h2r,h2i;
      i = 1;
      for (j=1;j<=n;j++)  {
        fft1[i++] = data1[j];
        fft1[i++] = data2[j];
      }
      FFT ( fft1,n,true );
      fft2[1] = fft1[2];    fft2[2] = 0.0;
      fft1[2] = 0.0;
      n2 = n + 2;
      for (j=2;j<=n/2+1;j++)  {
        bj = 2*j-1;    bn = 2*(n2-j)-1;
        h1r = 0.5*(fft1[bj]   + fft1[bn]);
        h1i = 0.5*(fft1[bj+1] - fft1[bn+1]);
        h2r = 0.5*(fft1[bj+1] + fft1[bn+1]);
        h2i = 0.5*(fft1[bn]   - fft1[bj]);
        fft1[bj] = h1r;    fft1[bj+1] = h1i;
        fft1[bn] = h1r;    fft1[bn+1] = -h1i;
        fft2[bj] = h2r;    fft2[bj+1] = h2i;
        fft2[bn] = h2r;    fft2[bn+1] = -h2i;
      }
    }

    void Convolve ( rvector data, int n, rvector respns, int m,
                    rvector ans,  bool Conv )  {
    //  Convolves or Deconvolves a real data set data[1:n] (including
    // any user-supplied zero padding) with a response function
    // respns[1..n], stored in wrap-around order in a real array of
    // length m<n (m should be an odd (3,5,7...) integer).  Wrap-around
    // order means that the first half of the array contains the impulse
    // response function at positive times, while the second half of
    // the array contains the impulse response function at negative
    // times, counting down from the highest element respns[m]. On
    // input Conv=true for convolution, false for deconvolution.
    // The answer is returned in the first n component of ans.
    // However, ans must be supplied in the calling program with
    // length at least 2*n, for consistency with TwoFFT.  n MUST
    // be an integer power of 2.
    //
    int      i,no2,rp,ip;
    rvector  fft;
    realtype B,D;

      GetVectorMemory ( fft,2*n,1 );
      for (i=1;i<=(m-1)/2;i++)
        respns[n+1-i] = respns[m+1-i];
      for (i=(m+3)/2;i<=n-(m-1)/2;i++)
        respns[i] = 0.0;

      TwoFFT ( data,respns,fft,ans,n );

      no2 = n/2;
      rp  = 1;   // pointer to real part
      ip  = 2;   // pointer to imaginary part
      for (i=1;i<=no2+1;i++)  {
        if (Conv) {
          B       = (fft[rp]*ans[rp] - fft[ip]*ans[ip])/no2;
          ans[ip] = (fft[ip]*ans[rp] + fft[rp]*ans[ip])/no2;
          ans[rp] = B;
        } else  {
          D = (ans[rp]*ans[rp] + ans[ip]*ans[ip])*no2;
          if (D==0.0)  {
            //  poor deconvolve at zero response
            ans[rp] = 0.0;
            ans[ip] = 0.0;
          } else  {
            B       = (fft[rp]*ans[rp] + fft[ip]*ans[ip])/D;
            ans[ip] = (fft[ip]*ans[rp] - fft[rp]*ans[ip])/D;
            ans[rp] = B;
          }
        }
        rp += 2;
        ip += 2;
      }

      ans[2] = ans[2*no2+1];
      FreeVectorMemory ( fft,1 );

      RealFFT ( ans,n,false );

    }


    void mConvolve ( rvector data, int n, int m )  {
    //
    //   Replaces array data[0..n-1] with the result of m recursive
    // convolutions (m>1) defined as
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
    realtype R,G,phi,B,m2,n2,d1;
    int      i,m1;

      if (m<1)  return;

      RealFFT ( data-1,n,true );

      m1 = m+1;
      m2 = m1/2.0;
      n2 = 2.0/n;
      d1 = data[1];
      for (i=0;i<=n;i+=2)  {
        if (i<n)  {
          R = data[i];
          if (i>1) G = data[i+1];
              else G = 0.0;
        } else  {
          R = d1;
          G = 0.0;
        }
        phi = atan2(G,R) * m1;
        B   = pow(R*R+G*G,m2);
        R   = B*cos(phi);
        G   = B*sin(phi);
        if (i<n)  {
          data[i]   = R*n2;
          data[i+1] = G*n2;
        } else
          data[1] = R*n2;
      }

      RealFFT ( data-1,n,false );

    }

  }  // namespace math

}  // namespace mmdb
