//  $Id: mmdb_math_graph.h $
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
//  **** Module  :  mmdb_math_graph  <interface>
//       ~~~~~~~~~
//  **** Namespace: mmdb::math::
//       ~~~~~~~~~~
//  **** Classes :  Vertex     ( graph vertex                        )
//       ~~~~~~~~~  Edge       ( graph edge                          )
//                  Graph      ( structural graph                    )
//                  Match      ( match of structural graphs          )
//                  GraphMatch ( CSIA algorithms for graphs matching )
//
//   (C) E. Krissinel 2000-2013
//
//  When used, please cite:
//
//   Krissinel, E. and Henrick, K. (2004)
//   Common subgraph isomorphism detection by backtracking search.
//   Software - Practice and Experience, 34, 591-607.
//
//  =================================================================
//

#ifndef  __MMDB_MATH_Graph__
#define  __MMDB_MATH_Graph__

#include <time.h>

#include "mmdb_atom.h"

namespace mmdb  {

  namespace math  {

    //  =========================  Vertex  ==========================

    DefineClass(Vertex);

    enum GRAPH_FLAG  {
      CHIRAL_RIGHT  = 0x10000000,
      CHIRAL_LEFT   = 0x20000000,
      ATOM_LEAVING  = 0x40000000,
      HYDROGEN_BOND = 0x0F000000,
      SYMREL_MASK   = 0x00FF0000,
      CHIRAL_MASK   = 0xCFFFFFFF,
      TYPE_MASK     = 0x00FFFFFF
    };

    class Vertex : public io::Stream  {

      friend class Graph;
      friend class GraphMatch;

      public:

        Vertex ();
        Vertex ( io::RPStream Object );
        Vertex ( int  vtype, cpstr vname );
        Vertex ( int  vtype );
        Vertex ( cpstr chem_elem );
        Vertex ( cpstr chem_elem, cpstr name );
        ~Vertex();

        void  SetVertex  ( cpstr chem_elem );
        void  SetVertex  ( int vtype, cpstr vname );
        void  SetVertex  ( int vtype   );
        void  SetType    ( int vtype   );
        void  SetTypeExt ( int typeExt );

        void  RemoveChirality();
        void  LeaveChirality ( int eltype );

        void  SetName     ( cpstr vname );
        void  SetProperty ( int vprop  );
        void  SetID       ( int vid    );
        void  AddBond     ();
        void  CopyNBonds  ( PVertex V );
        inline void  SetUserID   ( int vid ) { user_id = vid; }
        inline int   GetProperty () { return property; }
        inline int   GetID       () { return id;       }
        inline int   GetUserID   () { return user_id;  }
        inline cpstr GetName     () { return name;     }
        inline int   GetType     () { return type;     }
        inline int   GetTypeExt  () { return type_ext; }
        int   GetNBonds   ();

        void  SaveType    ();  // in userid
        void  RestoreType ();  // from userid
        void  CopyType    ( PVertex V );

        virtual void Print ( int PKey );

        virtual void Copy ( PVertex V );

        void  read  ( io::RFile f );
        void  write ( io::RFile f );

        void  mem_read  ( cpstr S, int & l );
        void  mem_write ( pstr S, int & l );

      protected:
        pstr name;     // name may be general, "C", "Hg", "Cl" etc.
        int  type;     // type of vertex, see comments in
                       // mmdb_math_graph.cpp
        int  type_ext; // vertex type extention
        int  property; // flagwise properties -- user-defined
        int  id;       // a graph-defined vertex id
        int  user_id;  // a user-defined vertex id

        void InitVertex();

    };

    DefineStreamFunctions(Vertex);


    //  ==========================  Edge  ===========================

    enum GRAPH_BOND  {
      BOND_SINGLE   = 1,
      BOND_DOUBLE   = 2,
      BOND_AROMATIC = 3,
      BOND_TRIPLE   = 4
    };

    DefineClass(Edge);

    class Edge : public io::Stream  {

      friend class Graph;
      friend class CGMatch;

      public:

        Edge ();
        Edge ( io::RPStream Object );
        Edge ( int vx1, int vx2, int btype );  // vx1,vx2 are numbered
                                               // as 1,2,3 on and refer
                                               // to vertices in the order
                                               // as they were added to
                                               // the graph; btype>0
        ~Edge();

        void  SetEdge ( int vx1, int vx2, cpstr btype );
        void  SetEdge ( int vx1, int vx2, int  btype ); // btype>0

        void  SetType     ( int btype );
        void  SetProperty ( int eprop );
        void  SaveType    ();  // in property
        void  RestoreType ();  // from property

        inline int GetVertex1  () { return v1;       }
        inline int GetVertex2  () { return v2;       }
        inline int GetType     () { return type;     }
        inline int GetProperty () { return property; }

        virtual void Print ( int PKey );

        virtual void Copy  ( PEdge G );

        void  read  ( io::RFile f );
        void  write ( io::RFile f );

        void  mem_read  ( cpstr S, int & l );
        void  mem_write ( pstr  S, int & l );

      protected:
        int  v1,v2;  //  >=1
        int  type;
        int  property;

        void  InitEdge();

    };

    DefineStreamFunctions(Edge);


    //  ==========================  Graph  ============================

    enum GRAPH_RC  {
      MKGRAPH_Ok            =  0,
      MKGRAPH_NoAtoms       = -1,
      MKGRAPH_ChangedAltLoc =  1,
      MKGRAPH_MaxOccupancy  =  2
    };

    DefineClass(Graph);

    class Graph : public io::Stream  {

      friend class GraphMatch;
      friend class CSBase0;

      public :

        Graph ();
        Graph ( PResidue R, cpstr altLoc=NULL );
        Graph ( io::RPStream Object );
        ~Graph();

        void  Reset   ();
        void  SetName ( cpstr gname );
        inline pstr GetName() { return name; }

        //   AddVertex(..) and AddEdge(..) do not copy the objects, but
        // take them over. This means that application should forget
        // about pointers to V and G once they were given to Graph.
        // Vertices and edges  must be allocated newly prior each call
        // to AddVertex(..) and AddEdge(..).
        void  AddVertex   ( PVertex  V );
        void  AddEdge     ( PEdge    G );
        void  SetVertices ( PPVertex V, int vlen );
        void  SetEdges    ( PPEdge   G, int glen );

        void  RemoveChirality();
        void  LeaveChirality ( int eltype );

        //   MakeGraph(..) makes a graph corresponding to residue R.
        // The graphs vertices then correspond to the residue's atoms
        // (Vertex::userid points to atom R->atom[Vertex::userid]),
        // edges are calculated as chemical bonds between atoms basing
        // on the table of cut-off distances.
        //   altCode specifies a particular conformation that should be
        // used for making the graph. If it is set to "" or NULL ("empty"
        // altcode) but the residue does not have conformation which
        // contains *only* ""-altcode atoms, a conformation corresponding
        // to maximal occupancy will be used. The same will happen if
        // altcode information in residue is not correct, whatever altCode
        // is specified.
        //   After making the graph, Build(..) should be called as usual
        // before graph matching.
        //   Non-negative return means that graph has been made.
        // MakeGraph(..) may return:
        //   MKGRAPH_Ok             everything is Ok
        //   MKGRAPH_NoAtoms        residue does not have atoms, graph
        //                          is not made
        //   MKGRAPH_ChangedAltLoc  a different altcode was used because
        //                          the residue has only one altcode and
        //                          that is different of
        //   MKGRAPH_MaxOccupancy   a maximal-occupancy conformation has
        //                          been chosen because of default
        //                          ""-altcode supplied or incorrect
        //                          altcode information in the residue
        int   MakeGraph   ( PResidue R, cpstr altLoc=NULL );

        int   MakeGraph   ( PPAtom atom, int nAtoms );

        void  HideType    ( int bond_vx_type );
        void  ExcludeType ( int type );

        void  MakeSymmetryRelief ( bool noCO2 );
        void  IdentifyRings      ();
        int   IdentifyConnectedComponents();  // returns their number >= 1

        int   Build       ( bool bondOrder );  // returns 0 if Ok

        void  MakeVertexIDs      ();  // simply numbers vertices as 1.. on
        int   GetVertexID        ( int vertexNo );
        int   GetVertexNo        ( cpstr vname  );
        // GetBondedVertexID(..) works after MoveType(..)
        int   GetNBondedVertices ( int vertexNo );
        int   GetBondedVertexID  ( int vertexNo, int bond_vx_type,
                                   int bondNo );

        PVertex   GetVertex ( int vertexNo );  // 1<=vertexNo<=nVertices
        inline int GetNofVertices() { return nVertices; }

        PEdge    GetEdge    ( int edgeNo );    // 1<=edgeNo<=nEdges
        inline int GetNofEdges() { return nEdges;    }

        void  GetVertices ( PPVertex & V, int & nV );
        void  GetEdges    ( PPEdge   & E, int & nE );

        virtual void Print();
        void  Print1();

        virtual void Copy ( PGraph G );

        void  read  ( io::RFile f );
        void  write ( io::RFile f );

        void  mem_read  ( cpstr S, int & l );
        void  mem_write ( pstr S, int & l );

      protected :
        pstr     name;
        int      nVertices,nEdges, nAllVertices,nAllEdges;
        PPVertex vertex;
        PPEdge   edge;
        imatrix  graph;

        void  InitGraph ();
        void  FreeMemory();

        void  markConnected ( int vno, int cno );

      private :
        int  nVAlloc,nEAlloc,nGAlloc;

    };

    DefineStreamFunctions(Graph);


    //  =========================  GMatch  ==========================

    DefineClass(GMatch);
    DefineStreamFunctions(GMatch);

    class GMatch : public io::Stream  {

      friend class GraphMatch;

      public :

        GMatch ();
        GMatch ( io::RPStream Object );
        GMatch ( ivector FV1, ivector FV2, int nv, int n, int m );
        ~GMatch();

        // FV1[] and FV2[] are copied into internal buffers
        void SetMatch ( ivector FV1, ivector FV2, int nv, int n, int m );

        bool isMatch       ( ivector FV1, ivector FV2, int nv );
        bool isCombination ( ivector FV1, ivector FV2, int nv );

        // do not allocate or dispose FV1 and FV2 in application!
        void GetMatch ( ivector & FV1, ivector & FV2, int & nv,
                        realtype & p1, realtype & p2 );

        void read  ( io::RFile f );
        void write ( io::RFile f );

        void mem_read  ( cpstr S, int & l );
        void mem_write ( pstr  S, int & l );

      protected :
        int     n1,n2,mlength;
        ivector F1,F2;

        void InitGMatch();

      private :
        int nAlloc;

    };


    //  =======================  GraphMatch  =========================

    #define  _UseRecursion

    enum GRAPH_MATCH_FLAG  {
      GMF_UniqueMatch    = 0x00000001,
      GMF_NoCombinations = 0x00000002
    };

    enum VERTEX_EXT_TYPE  {
      EXTTYPE_Ignore   = 0,
      EXTTYPE_Equal    = 1,
      EXTTYPE_AND      = 2,
      EXTTYPE_OR       = 3,
      EXTTYPE_XOR      = 4,
      EXTTYPE_NotEqual = 5,
      EXTTYPE_NotAND   = 6,
      EXTTYPE_NotOR    = 7
    };

    DefineClass(GraphMatch);

    class GraphMatch : public io::Stream  {

      public :

        GraphMatch ();
        GraphMatch ( io::RPStream Object );
        ~GraphMatch();

        void SetFlag          ( word flag );
        void RemoveFlag       ( word flag );
        void SetMaxNofMatches ( int maxNofMatches, bool stopOnMaxN );
        void SetTimeLimit     ( int maxTimeToRun=0 );
        inline bool GetStopSignal() { return Stop; }

        void Reset();

        //  MatchGraphs looks for maximal common subgraphs of size
        //  not less than minMatch. The number of found subgraphs
        //  is returned by GetNofMatches(), the subgraph vertices
        //  are returned by GetMatch(..). Control parameters:
        //      vertexType   true if vertex type should be taken
        //                   into account and False otherwise
        //      vertexExt    key to use extended vertex types (defined
        //                   as type_ext in Vertex).
        void MatchGraphs    ( PGraph Gh1, PGraph Gh2, int minMatch,
                              bool vertexType=true,
                              VERTEX_EXT_TYPE vertexExt=EXTTYPE_Ignore );
        void PrintMatches   ();
        inline int GetNofMatches  () { return nMatches; }
        inline int GetMaxMatchSize() { return maxMatch; }

        // do not allocate or dispose FV1 and FV2 in application!
        // FV1/p1 will always correspond to Gh1, and FV2/p2 -
        // to Gh2 as specified in MatchGraphs(..)
        void GetMatch ( int MatchNo, ivector & FV1, ivector & FV2,
                        int & nv, realtype & p1, realtype & p2 );

        void read  ( io::RFile f );
        void write ( io::RFile f );

        void mem_read  ( cpstr S, int & l );
        void mem_write ( pstr  S, int & l );

      protected :
        PGraph   G1,G2;
        PPVertex V1;
        PPVertex V2;
        imatrix  c1,c2;
        bool     swap;
    #ifndef _UseRecursion
        ivector  jj;
    #endif
        int      n,m;

        imatrix3 P;
        imatrix  iF1;
        ivector  F1,F2,ix;

        int      nMatches,maxNMatches;
        PPGMatch Match;
        bool     wasFullMatch,Stop,stopOnMaxNMathches;
        word     flags;
        int      maxMatch,timeLimit;

        void  InitGraphMatch();
        void  FreeMemory    ();
        void  FreeRecHeap   ();
        void  GetMemory     ();
        void  GetRecHeap    ();
        int   Initialize    ( bool vertexType, int vertexExt );
    #ifdef _UseRecursion
        void  Backtrack     ( int i );          // exact matching
    #else
        void  Ullman        ();
    #endif
        void  Backtrack1    ( int i, int k0 );  // exact/partial matching
        void  CollectMatch  ( int nm );

      private :
        int     nAlloc,mAlloc,nMAlloc;
        time_t  startTime;

    };

    DefineStreamFunctions(GraphMatch);

    extern void  SetGraphAllocPortion ( int alloc_portion );

    /*
    extern void  TestGraphMatch();
    */

  }  // namespace math

}  // namespace mmdb

#endif
