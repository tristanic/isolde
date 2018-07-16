//  $Id: mmdb_uddata.cpp $
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
//  **** Module  :   MMDBF_UDData <implementation>
//       ~~~~~~~~~
//  **** Project :   MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//
//  **** Classes :   mmdb::UDData ( user-defined data )
//       ~~~~~~~~~
//
//   (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#include <string.h>

#include "mmdb_uddata.h"

namespace mmdb  {

  //  ========================  UDRegister  ==========================

  #define  nUDRTypes  5

  UDRegister::UDRegister() : io::Stream()  {
    InitUDRegister();
  }

  UDRegister::UDRegister ( io::RPStream Object ) : io::Stream(Object)  {
    InitUDRegister();
  }

  UDRegister::~UDRegister()  {
    FreeUDRegister();
  }

  void  UDRegister::InitUDRegister()  {
  int i;
    for (i=0;i<nUDRTypes;i++)  {
      nIUDR[i]       = 0;
      nRUDR[i]       = 0;
      nSUDR[i]       = 0;
      IUDRegister[i] = NULL;
      RUDRegister[i] = NULL;
      SUDRegister[i] = NULL;
    }
  }

  void  UDRegister::FreeUDRegister()  {
  int i,j;

    for (j=0;j<nUDRTypes;j++)  {

      if (IUDRegister[j])  {
        for (i=0;i<nIUDR[j];i++)
          if (IUDRegister[j][i])  delete[] IUDRegister[j][i];
        delete[] IUDRegister[j];
        IUDRegister[j] = NULL;
      }
      nIUDR[j] = 0;

      if (RUDRegister[j])  {
        for (i=0;i<nRUDR[j];i++)
          if (RUDRegister[j][i])  delete[] RUDRegister[j][i];
        delete[] RUDRegister[j];
        RUDRegister[j] = NULL;
      }
      nRUDR[j] = 0;
      if (SUDRegister[j])  {
        for (i=0;i<nRUDR[j];i++)
          if (SUDRegister[j][i])  delete[] SUDRegister[j][i];
        delete[] SUDRegister[j];
        SUDRegister[j] = NULL;
      }
      nSUDR[j] = 0;

    }

  }

  int UDRegister::RegisterUDData ( psvector & UDRegister,
                                   int      & nUDR,
                                   cpstr      UDDataID )  {
  psvector UDReg;
  int      i,UDDhandle,n;

    n         = -1;
    UDDhandle = 0;
    for (i=0;(i<nUDR) && (!UDDhandle);i++)
      if (UDRegister[i])  {
        if (!strcmp(UDDataID,UDRegister[i]))
          UDDhandle = i+1;
      } else
        n = i;

    if (!UDDhandle)  {
      if (n<0)  {
        UDReg = new pstr[nUDR+1];
        for (i=0;i<nUDR;i++)
          UDReg[i] = UDRegister[i];
        UDReg[nUDR] = NULL;
        if (UDRegister)  delete[] UDRegister;
        UDRegister = UDReg;
        n = nUDR;
        nUDR++;
      }
      CreateCopy ( UDRegister[n],UDDataID );
      UDDhandle = n+1;
    }

    return UDDhandle;

  }

  static int UDRegisterFlag[nUDRTypes] = {
    UDRF_ATOM,
    UDRF_RESIDUE,
    UDRF_CHAIN,
    UDRF_MODEL,
    UDRF_HIERARCHY
  };


  int  UDRegister::RegisterUDInteger ( UDR_TYPE udr_type,
                                       cpstr UDDataID )  {
    if ((udr_type>=0) && (udr_type<nUDRTypes))
      return RegisterUDData ( IUDRegister[udr_type],
                              nIUDR[udr_type],UDDataID ) |
             UDRegisterFlag[udr_type];
    else
      return UDDATA_WrongUDRType;
  }

  int  UDRegister::RegisterUDReal ( UDR_TYPE udr_type,
                                    cpstr UDDataID )  {
    if ((udr_type>=0) && (udr_type<nUDRTypes))
      return RegisterUDData ( RUDRegister[udr_type],
                              nRUDR[udr_type],UDDataID ) |
             UDRegisterFlag[udr_type];
    else
      return UDDATA_WrongUDRType;
  }

  int  UDRegister::RegisterUDString ( UDR_TYPE udr_type,
                                      cpstr UDDataID )  {
    if ((udr_type>=0) && (udr_type<nUDRTypes))
      return RegisterUDData ( SUDRegister[udr_type],
                              nSUDR[udr_type],UDDataID ) |
             UDRegisterFlag[udr_type];
    else
      return UDDATA_WrongUDRType;
  }

  int  UDRegister::GetUDDHandle ( UDR_TYPE udr_type,
                                  cpstr UDDataID )  {
  int  i,UDDhandle;

    if ((udr_type>=0) && (udr_type<nUDRTypes))  {

      UDDhandle = 0;

      for (i=0;(i<nIUDR[udr_type]) && (!UDDhandle);i++)
        if (IUDRegister[udr_type][i])  {
          if (!strcmp(UDDataID,IUDRegister[udr_type][i]))
            UDDhandle = i+1;
        }
      for (i=0;(i<nRUDR[udr_type]) && (!UDDhandle);i++)
        if (RUDRegister[udr_type][i])  {
          if (!strcmp(UDDataID,RUDRegister[udr_type][i]))
            UDDhandle = i+1;
        }
      for (i=0;(i<nSUDR[udr_type]) && (!UDDhandle);i++)
        if (SUDRegister[udr_type][i])  {
          if (!strcmp(UDDataID,SUDRegister[udr_type][i]))
            UDDhandle = i+1;
        }

      if (UDDhandle)  return UDDhandle | UDRegisterFlag[udr_type];
                else  return UDDhandle;

    } else
      return UDDATA_WrongUDRType;

  }


  void  UDRegister::write ( io::RFile f )  {
  int  i,j;
  byte Version=1;
    f.WriteByte ( &Version );
    for (j=0;j<nUDRTypes;j++)  {
      f.WriteInt ( &nIUDR[j] );
      for (i=0;i<nIUDR[j];i++)
        f.CreateWrite ( IUDRegister[j][i] );
      f.WriteInt ( &nRUDR[j] );
      for (i=0;i<nRUDR[j];i++)
        f.CreateWrite ( RUDRegister[j][i] );
      f.WriteInt ( &nSUDR[j] );
      for (i=0;i<nSUDR[j];i++)
        f.CreateWrite ( SUDRegister[j][i] );
    }
  }

  void  UDRegister::read ( io::RFile f )  {
  int  i,j;
  byte Version;
    f.ReadByte ( &Version );
    FreeUDRegister();
    for (j=0;j<nUDRTypes;j++)  {
      f.ReadInt ( &nIUDR[j] );
      if (nIUDR[j]>0)  {
        IUDRegister[j] = new pstr[nIUDR[j]];
        for (i=0;i<nIUDR[j];i++)  {
          IUDRegister[j][i] = NULL;
          f.CreateRead ( IUDRegister[j][i] );
        }
      }
      f.ReadInt ( &nRUDR[j] );
      if (nRUDR[j]>0)  {
        RUDRegister[j] = new pstr[nRUDR[j]];
        for (i=0;i<nRUDR[j];i++)  {
          RUDRegister[j][i] = NULL;
          f.CreateRead ( RUDRegister[j][i] );
        }
      }
      f.ReadInt ( &nSUDR[j] );
      if (nSUDR[j]>0)  {
        SUDRegister[j] = new pstr[nSUDR[j]];
        for (i=0;i<nSUDR[j];i++)  {
          SUDRegister[j][i] = NULL;
          f.CreateRead ( SUDRegister[j][i] );
        }
      }
    }
  }


  MakeStreamFunctions(UDRegister)



  //  ==========================  UDData  ============================

  UDData::UDData() : Mask()  {
    InitUDData();
  }

  UDData::UDData ( io::RPStream Object ) : Mask(Object)  {
    InitUDData();
  }

  UDData::~UDData()  {
    FreeUDDMemory();
  }

  void UDData::InitUDData()  {
    IUData = NULL;
    RUData = NULL;
    SUData = NULL;
  }

  void UDData::FreeUDDMemory()  {
  int i,l;
    FreeVectorMemory ( IUData,0 );
    FreeVectorMemory ( RUData,0 );
    if (SUData)  {
      l = getNofSUData();
      for (i=0;i<=l;i++)
        if (SUData[i])  delete[] SUData[i];
      delete[] SUData;
    }
    IUData = NULL;
    RUData = NULL;
    SUData = NULL;
  }

  int  UDData::getNofIUData()  {
    if (!IUData)  return 0;
    return IUData[0];
  }

  int  UDData::getNofRUData()  {
    if (!RUData)  return 0;
    return mround(RUData[0]);
  }

  int  UDData::getNofSUData()  {
    if (!SUData)    return 0;
    if (!SUData[0]) return 0;
    return (int(SUData[0][0]) << 24) +
           (int(SUData[0][1]) << 16) +
           (int(SUData[0][2]) << 8)  +
            int(SUData[0][3]);
  }

  void UDData::setNofSUData ( int newN )  {
    if (!SUData)    return;
    if (!SUData[0]) return;
    SUData[0][3] = byte( newN & 0x000000FF);
    SUData[0][2] = byte((newN & 0x0000FF00) >> 8);
    SUData[0][1] = byte((newN & 0x00FF0000) >> 16);
    SUData[0][0] = byte((newN & 0xFF000000) >> 24);
  }


  int  UDData::putUDData ( int UDDhandle, int iudd )  {
  ivector IUD;
  int     i,l,udh;
    udh = UDDhandle & UDRF_MASK;
    if (udh<1)  return UDDATA_WrongHandle;
    l = getNofIUData();
    if (udh>l)  {
      GetVectorMemory ( IUD,udh+1,0 );
      IUD[0] = udh;
      for (i=1;i<=l;i++)
        IUD[i] = IUData[i];
      for (i=l+1;i<udh;i++)
        IUD[i] = MinInt4;
      FreeVectorMemory ( IUData,0 );
      IUData = IUD;
    }
    IUData[udh] = iudd;
    return UDDATA_Ok;
  }

  int  UDData::putUDData ( int UDDhandle, realtype rudd )  {
  rvector RUD;
  int     i,l,udh;
    udh = UDDhandle & UDRF_MASK;
    if (udh<1)  return UDDATA_WrongHandle;
    l = getNofRUData();
    if (udh>l)  {
      GetVectorMemory ( RUD,udh+1,0 );
      RUD[0] = udh;
      for (i=1;i<=l;i++)
        RUD[i] = RUData[i];
      for (i=l+1;i<udh;i++)
        RUD[i] = -MaxReal;
      FreeVectorMemory ( RUData,0 );
      RUData = RUD;
    }
    RUData[udh] = rudd;
    return UDDATA_Ok;
  }

  int  UDData::putUDData ( int UDDhandle, cpstr sudd )  {
  psvector SUD;
  int      i,l,udh;
    udh = UDDhandle & UDRF_MASK;
    if (udh<1)  return UDDATA_WrongHandle;
    l = getNofSUData();
    if (udh>l)  {
      if (l>0)  {
        GetVectorMemory ( SUD,udh+1,0 );
        for (i=0;i<=l;i++)
          SUD[i] = SUData[i];
        for (i=l+1;i<=udh;i++)
          SUD[i] = NULL;
        FreeVectorMemory ( SUData,0 );
        SUData = SUD;
      } else  {
        GetVectorMemory ( SUData,udh+1,0 );
        SUData[0] = new char[4];
        for (i=1;i<=udh;i++)
          SUData[i] = NULL;
      }
      setNofSUData ( udh );
    }
    CreateCopy ( SUData[udh],sudd );
    return UDDATA_Ok;
  }

  int  UDData::getUDData ( int UDDhandle, int & iudd )  {
  int l,udh;
    iudd = 0;
    udh  = UDDhandle & UDRF_MASK;
    if (udh<1)  return UDDATA_WrongHandle;
    l = getNofIUData();
    if (udh>l)  return UDDATA_NoData;
    iudd = IUData[udh];
    if (iudd==MinInt4)  return UDDATA_NoData;
    return UDDATA_Ok;
  }

  int  UDData::getUDData ( int UDDhandle, realtype & rudd )  {
  int l,udh;
    rudd = 0.0;
    udh = UDDhandle & UDRF_MASK;
    if (udh<1)  return UDDATA_WrongHandle;
    l = getNofRUData();
    if (udh>l)  return UDDATA_NoData;
    rudd = RUData[udh];
    if (rudd==-MaxReal)  return UDDATA_NoData;
    return UDDATA_Ok;
  }

  int  UDData::getUDData ( int UDDhandle, pstr sudd, int maxLen )  {
  int l,udh;
    sudd[0] = char(0);
    udh = UDDhandle & UDRF_MASK;
    if (udh<1)  return UDDATA_WrongHandle;
    l = getNofSUData();
    if (udh>l)  return UDDATA_NoData;
    if (!SUData[udh])  return UDDATA_NoData;
    strcpy_n0 ( sudd,SUData[udh],maxLen-1 );
    return UDDATA_Ok;
  }

  pstr  UDData::getUDData ( int UDDhandle, int * retcode )  {
  int l,udh;
    udh = UDDhandle & UDRF_MASK;
    if (udh<1)  {
      if (retcode)  *retcode = UDDATA_WrongHandle;
      return NULL;
    }
    l = getNofSUData();
    if (udh>l)  {
      if (retcode)  *retcode = UDDATA_NoData;
      return NULL;
    }
    if (!SUData[udh])  {
      if (retcode)  *retcode = UDDATA_NoData;
      return NULL;
    }
    if (retcode)  *retcode = UDDATA_Ok;
    return SUData[udh];
  }

  int  UDData::getUDData ( int UDDhandle, pstr & sudd )  {
  int l,udh;
    udh = UDDhandle & UDRF_MASK;
    if (udh<1)  {
      if (sudd)  {
        delete[] sudd;
        sudd = NULL;
      }
      return UDDATA_WrongHandle;
    }
    l = getNofSUData();
    if (udh>l)  {
      if (sudd)  {
        delete[] sudd;
        sudd = NULL;
      }
      return UDDATA_NoData;
    }
    if (!SUData[udh])  {
      if (sudd)  {
        delete[] sudd;
        sudd = NULL;
      }
      return UDDATA_NoData;
    }
    CreateCopy ( sudd,SUData[udh] );
    return UDDATA_Ok;
  }


  void  UDData::write ( io::RFile f )  {
  int  i,l;
  byte Version=1;

    f.WriteByte ( &Version );

    Mask::write ( f );

    if (IUData)  l = IUData[0];
           else  l = -1;
    f.WriteVector ( IUData,l+1,0 );
    if (RUData)  l = mround(RUData[0]);
           else  l = -1;
    f.WriteVector ( RUData,l+1,0 );
    l = getNofSUData();
    f.WriteInt ( &l );
    for (i=1;i<=l;i++)
      f.CreateWrite ( SUData[i] );
  }

  void  UDData::read  ( io::RFile f )  {
  int  i,l;
  byte Version;

    f.ReadByte ( &Version );

    FreeUDDMemory();

    Mask::read ( f );

    f.CreateReadVector ( IUData,0 );
    f.CreateReadVector ( RUData,0 );
    f.ReadInt ( &l );
    if (l>0)  {
      SUData = new pstr[l+1];
      SUData[0] = new char[4];
      setNofSUData ( l );
      for (i=1;i<=l;i++)  {
        SUData[i] = NULL;
        f.CreateRead ( SUData[i] );
      }
    }
  }


  MakeStreamFunctions(UDData)

}  // namespace mmdb

