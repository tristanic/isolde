//  $Id: mmdb_bondmngr.cpp $
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
//    15.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  mmdb_bondmngr <implementation>
//       ~~~~~~~~~
//       Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::BondManager ( MMDB bonds maker )
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//


#include <string.h>

#include "mmdb_bondmngr.h"
#include "mmdb_math_graph.h"

namespace mmdb  {

  //  =====================   BondManager   =====================

  BondManager::BondManager() : SelManager()  {
  }

  BondManager::BondManager ( io::RPStream Object )
             : SelManager(Object)  {
  }

  BondManager::~BondManager()  {}

  void  BondManager::MakeBonds ( bool calc_only )  {
  UNUSED_ARGUMENT(calc_only);
  PModel         mdl;
  PChain         chain;
  PResidue       res;
  math::Graph    graph;
  math::PPVertex V;
  math::PPEdge   E;
  int            i, im,ic,ir, nV,nE, k1,k2;

    RemoveBonds();

    for (im=0;im<nModels;im++)  {
      mdl = model[im];
      if (mdl)
        for (ic=0;ic<mdl->nChains;ic++)  {
          chain = mdl->chain[ic];
          if (chain)
            for (ir=0;ir<chain->nResidues;ir++)  {
              res = chain->residue[ir];
              if (res)  {
                graph.MakeGraph   ( res,NULL );
                graph.GetVertices ( V,nV );
                graph.GetEdges    ( E,nE );
                for (i=0;i<nE;i++)  {
                  k1 = V[E[i]->GetVertex1()-1]->GetUserID();
                  k2 = V[E[i]->GetVertex2()-1]->GetUserID();
                  res->atom[k1]->AddBond ( res->atom[k2],E[i]->GetType() );
                  res->atom[k2]->AddBond ( res->atom[k1],E[i]->GetType() );
                }
              }
            }
        }
    }

  }

  void  BondManager::RemoveBonds()  {
  int i;
    for (i=0;i<nAtoms;i++)
      if (atom[i])
        atom[i]->FreeBonds();
  }

  //  -------------------  Stream functions  ----------------------

  void  BondManager::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte ( &Version );
    SelManager::write ( f );
  }

  void  BondManager::read ( io::RFile f )  {
  byte Version;
    f.ReadByte ( &Version );
    SelManager::read ( f );
  }


  MakeStreamFunctions(BondManager)

}  // namespace mmdb
