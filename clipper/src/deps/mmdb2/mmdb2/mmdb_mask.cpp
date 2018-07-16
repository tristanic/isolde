//  $Id: mmdb_mask.cpp $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2013.
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
//  **** Module  :   MMDBF_Mask <implementation>
//       ~~~~~~~~~
//  **** Project :   MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//
//  **** Classes :   mmdb::Mask  ( atom selection mask )
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#include <string.h>
#include <stdlib.h>

#include "mmdb_mask.h"

namespace mmdb  {

  //  ====================  Mask  ========================

  Mask::Mask() : io::Stream()  {
    InitMask();
  }

  Mask::Mask ( io::RPStream Object ) : io::Stream(Object)  {
    InitMask();
  }

  Mask::~Mask()  {
    ClearMask();
  }

  void Mask::InitMask()  {
    mlen = 0;
    m    = NULL;
  }

  void Mask::SetMaskBit ( int BitNo )  {
  int n,i;
    n = BitNo/(8*sizeof(word));
    Expand ( n+1 );
    i = BitNo - n*(8*sizeof(word));
    m[n] |= ((word)1 << i);
  }

  void Mask::Expand ( int n )  {
  wvector m1;
  int     i;
    if (mlen<n)  {
      m1 = new word[n];
      for (i=0;i<mlen;i++)
        m1[i] = m[i];
      for (i=mlen;i<n;i++)
        m1[i] = 0;
      if (m)  delete[] m;
      m    = m1;
      mlen = n;
    }
  }

  void  Mask::NewMask ( PPMask Mask, int nMasks )  {
  int  i,nlen;
  word w;
    ClearMask();
    if (Mask && (nMasks>0))  {
      nlen = 0;
      w    = 0;
      while (w==0)  {
        for (i=0;i<nMasks;i++)
          if (Mask[i])  {
            if (nlen<Mask[i]->mlen)
              w |= Mask[i]->m[nlen];
          }
        nlen++;
        w = ~w;
      }
      Expand ( nlen );
      i    = nlen-1;
      m[i] = 1;
      while (!(m[i] & w))
        m[i] <<= 1;
    } else  {
      Expand ( 1 );
      m[0] = 1;
    }
  }

  void  Mask::CopyMask ( PMask Mask )  {
  int i;
    if (mlen!=Mask->mlen)  ClearMask();
    if (Mask)  {
      mlen = Mask->mlen;
      if (mlen>0)  {
        m = new word[mlen];
        for (i=0;i<mlen;i++)
          m[i] = Mask->m[i];
      }
    }
  }

  void  Mask::SetMask ( PMask Mask )  {
  int i;
    if (Mask) {
      Expand ( Mask->mlen );
      for (i=0;i<Mask->mlen;i++)
        m[i] |= Mask->m[i];
    }
  }

  void  Mask::RemoveMask ( PMask Mask )  {
  int i,l;
    if (Mask) {
      l = IMin(mlen,Mask->mlen);
      for (i=0;i<l;i++)
        m[i] &= ~Mask->m[i];
    }
  }

  void  Mask::SelMask ( PMask Mask )  {
  int i,l;
    if (Mask)  {
      l = IMin(mlen,Mask->mlen);
      for (i=0;i<l;i++)
        m[i] &= Mask->m[i];
      for (i=l;i<mlen;i++)
        m[i] = 0;
    } else
      ClearMask();
  }

  void  Mask::XadMask ( PMask Mask )  {
  int i;
    if (Mask) {
      Expand ( Mask->mlen );
      for (i=0;i<Mask->mlen;i++)
        m[i] ^= Mask->m[i];
    }
  }

  void  Mask::ClearMask()  {
    if (m)  delete[] m;
    m = NULL;
    mlen = 0;
  }

  void  Mask::NegMask()  {
  int i;
    for (i=0;i<mlen;i++)
      m[i] = ~m[i];
  }

  bool  Mask::CheckMask ( PMask Mask )  {
  int i,l;
    if (Mask)  {
      i = 0;
      l = IMin(mlen,Mask->mlen);
      while ((i<l) && (!(m[i] & Mask->m[i])))  i++;
      return (i<l);
    } else
      return false;
  }

  bool  Mask::isMask()  {
  int i=0;
    while ((i<mlen) && (!m[i]))  i++;
    return (i<mlen);
  }

  pstr  Mask::Print ( pstr S )  {
  int  i,j,k;
  word w;
    j = 0;
    for (i=0;i<mlen;i++)  {
      w = 1;
      for (k=0;k<8*(int)sizeof(word);k++)  {
        if (w & m[i])  S[j] = '1';
                 else  S[j] = '0';
        w <<= 1;
        j++;
      }
    }
    S[j] = char(0);
    return S;
  }

  void  Mask::write ( io::RFile f )  {
  int i;
    f.WriteInt ( &mlen );
    for (i=0;i<mlen;i++)
      f.WriteWord ( &(m[i]) );
  }

  void  Mask::read ( io::RFile f )  {
  int i;
    if (m)  {
      delete[] m;
      m = NULL;
    }
    f.ReadInt ( &mlen );
    if (mlen>0)  {
      m = new word[mlen];
      for (i=0;i<mlen;i++)
        f.ReadWord ( &(m[i]) );
    }
  }

  MakeStreamFunctions(Mask)

}  // namespace mmdb

