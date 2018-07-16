//  $Id: mmdb_xml_.h $
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
//    06.12.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  MMDB_XML <interface>
//       ~~~~~~~~~
//  **** Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::xml::XMLObject
//       ~~~~~~~~~
//
//   (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#ifndef __MMDB_XML__
#define __MMDB_XML__

#include "mmdb_mmcif_.h"

namespace mmdb  {

  namespace xml  {


    //  ======================  XMLObject  ==========================

    enum XML_RC  {
      XMLRC_Ok           = 0,
      XMLRC_NoFile       = 1,
      XMLRC_CantOpenFile = 2,
      XMLRC_NoTag        = 3,
      XMLRC_BrokenTag    = 4,
      XMLRC_UnclosedTag  = 5,
      XMLRC_RFormatError = 6,
      XMLRC_IFormatError = 7,
      XMLRC_OFormatError = 8
    };

    DefineClass(XMLObject);
    DefineStreamFunctions(XMLObject);

    class XMLObject : public io::Stream  {

      public :

        XMLObject ();
        XMLObject ( cpstr Tag );
        XMLObject ( cpstr Tag, cpstr Data );
        XMLObject ( cpstr Tag, realtype V, int length=11 );
        XMLObject ( cpstr Tag, int     iV, int length=0  );
        XMLObject ( cpstr Tag, bool    bV );
        XMLObject ( cpstr Tag, PXMLObject XMLObject );
        XMLObject ( io::RPStream Object );
        ~XMLObject();

        void  SetTag       ( cpstr Tag                             );
        void  AddAttribute ( cpstr name, cpstr      value          );
        void  AddAttribute ( cpstr name, const int  iV             );
        void  AddAttribute ( cpstr name, const bool bV             );
        void  SetData      ( cpstr Data                            );
        void  AddData      ( cpstr Data                            );
        void  SetData      ( const realtype V, const int length=11 );
        void  SetData      ( const int     iV, const int length=0  );
        void  SetData      ( const bool    bV                      );

        int AddMMCIFCategory ( mmcif::PCategory mmCIFCat    );
        int AddMMCIFStruct   ( mmcif::PStruct   mmCIFStruct );
        int AddMMCIFLoop     ( mmcif::PLoop     mmCIFLoop   );
        int AddMMCIFData     ( mmcif::PData     mmCIFData   );

        inline pstr GetTag()  { return objTag; }

        //   Here and below the functions allow for "tag1>tag2>tag3>..."
        // as a composite multi-level tag, e.g. the above may stand for
        // <tag1><tag2><tag3>data</tag3></tag2></tag1>. NULL tag
        // corresponds to "this" object.
        //   objNo counts same-tag objects of the *highest* level used
        // (e.g. level tag3 for composite tag  tag1>tag2>tag3 ).
        //   GetData ( pstr& ... ) only copies a pointer to data.
        pstr   GetData ( cpstr Tag=NULL, int objNo=1 );
        XML_RC GetData ( pstr   & Data, cpstr Tag=NULL, int objNo=1 );
        XML_RC GetData ( realtype &  V, cpstr Tag=NULL, int objNo=1 );
        XML_RC GetData ( int      & iV, cpstr Tag=NULL, int objNo=1 );
        XML_RC GetData ( bool     & bV, cpstr Tag=NULL, int objNo=1 );

        PXMLObject GetObject     ( cpstr Tag, int objNo=1 );
        PXMLObject GetFirstObject();
        PXMLObject GetLastObject ();
        inline int GetNumberOfObjects() { return nObjects; }
        PXMLObject GetObject     ( int objectNo );  // 0,1,...

        inline PXMLObject GetParent() { return parent; }

        void   AddObject    ( PXMLObject XMLObject, int lenInc=10 );
        void   InsertObject ( PXMLObject XMLObject, int pos,
                              int lenInc=10 );

        XML_RC WriteObject ( cpstr FName, int pos=0, int ident=2  );
        void   WriteObject ( io::RFile f, int pos=0, int ident=2  );
        XML_RC ReadObject  ( cpstr FName );
        XML_RC ReadObject  ( io::RFile f, pstr S, int & pos, int slen );

        virtual void Copy ( PXMLObject xmlObject );

        void  write ( io::RFile f );
        void  read  ( io::RFile f );

      protected:
        PXMLObject  parent;
        pstr        objTag;
        pstr        objData;
        int         nObjects,nAlloc;
        PPXMLObject object;
        int         nAttributes,nAttrAlloc;
        psvector    attr_name,attr_value;

        void         InitXMLObject();
        virtual void FreeMemory   ();

        inline void SetParent ( PXMLObject p ) { parent = p; }

    };


    extern PXMLObject mmCIF2XML ( mmcif::PData mmCIFData,
                                  int * rc=NULL );
    extern PXMLObject mmCIF2XML ( cpstr XMLName, mmcif::PFile mmCIFFile,
                                  int * rc=NULL );

  }  // namespace xml

}  // namespace mmdb

#endif


