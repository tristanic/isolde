//  $Id: mmdb_xml_.cpp $
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
//  **** Module  :  MMDB_XML <implementation>
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

#include <stdlib.h>
#include <string.h>

#include "mmdb_xml_.h"

namespace mmdb  {

  namespace xml  {

    //  ======================  XMLObject  ==========================

    XMLObject::XMLObject() : io::Stream()  {
      InitXMLObject();
    }

    XMLObject::XMLObject ( cpstr Tag ) : io::Stream()  {
      InitXMLObject();
      SetTag  ( Tag );
    }

    XMLObject::XMLObject ( cpstr Tag, cpstr Data ) : io::Stream()  {
      InitXMLObject();
      SetTag  ( Tag  );
      SetData ( Data );
    }

    XMLObject::XMLObject ( cpstr Tag, realtype V, int length )
              : io::Stream()  {
      InitXMLObject();
      SetTag  ( Tag      );
      SetData ( V,length );
    }

    XMLObject::XMLObject ( cpstr Tag, int iV, int length )
              : io::Stream()  {
      InitXMLObject();
      SetTag  ( Tag       );
      SetData ( iV,length );
    }

    XMLObject::XMLObject ( cpstr Tag, bool bV ) : io::Stream()  {
      InitXMLObject();
      SetTag  ( Tag );
      SetData ( bV  );
    }

    XMLObject::XMLObject ( cpstr Tag, PXMLObject XMLObject )
              : io::Stream()  {
      InitXMLObject();
      SetTag    ( Tag       );
      AddObject ( XMLObject );
    }


    XMLObject::XMLObject ( io::RPStream Object ) : io::Stream(Object)  {
      InitXMLObject();
    }

    XMLObject::~XMLObject()  {
      FreeMemory();
    }

    void  XMLObject::InitXMLObject()  {
      parent      = NULL;
      objTag      = NULL;
      objData     = NULL;
      nObjects    = 0;
      nAlloc      = 0;
      object      = NULL;
      nAttributes = 0;
      nAttrAlloc  = 0;
      attr_name   = NULL;
      attr_value  = NULL;
    }

    void  XMLObject::FreeMemory()  {
    int i;

      if (objTag)   delete[] objTag;
      if (objData)  delete[] objData;

      objTag   = NULL;
      objData  = NULL;

      if (object)  {
        for (i=0;i<nAlloc;i++)
          if (object[i])  delete object[i];
        delete[] object;
      }

      nObjects = 0;
      nAlloc   = 0;
      object   = NULL;

      if (attr_name)  {
        for (i=0;i<nAttrAlloc;i++)  {
          if (attr_name [i])  delete[] attr_name [i];
          if (attr_value[i])  delete[] attr_value[i];
        }
        FreeVectorMemory ( attr_name ,0 );
        FreeVectorMemory ( attr_value,0 );
      }

      nAttributes = 0;
      nAttrAlloc  = 0;
      attr_name   = NULL;
      attr_value  = NULL;

    }

    void  XMLObject::SetTag ( cpstr Tag )  {
    pstr p,t;
    int  n;

      // count ampersands
      p = pstr(Tag);
      n = 0;
      while (*p)  {
        if (*p=='&')  n++;
        p++;
      }
      // calculate the tag length
      n = n*4 + strlen(Tag) + 1;
      // substract leading underscores
      p = pstr(Tag);
      while (*p=='_')  {
        p++;
        n--;
      }
      // allocate tag space
      if (objTag)  delete[] objTag;
      objTag = new char[n];
      // copy tag, replacing square brackets and ampersands
      t = objTag;
      while (*p)  {
        if (*p=='[')  {
          *t = '-';
          t++;
        } else if (*p==']')  {
          if ((p[1]) && (p[1]!='['))  {
            *t = '-';
            t++;
          }
        } else if (*p=='&')  {
          strcpy ( t,"_and_" );
          if (p[1])  t += 5;
               else  t += 4;
        } else  {
          *t = *p;
          t++;
        }
        p++;
      }
      *t = char(0);
    }

    void  XMLObject::AddAttribute ( cpstr name, cpstr value )  {
    psvector an,av;
    int      i;

      if (nAttributes>=nAttrAlloc)  {
        nAttrAlloc = nAttributes + 10;
        GetVectorMemory ( an,nAttrAlloc,0 );
        GetVectorMemory ( av,nAttrAlloc,0 );
        for (i=0;i<nAttrAlloc;i++)  {
          an[i] = NULL;
          av[i] = NULL;
        }
        for (i=0;i<nAttributes;i++)  {
          CreateCopy ( an[i],attr_name [i] );
          CreateCopy ( av[i],attr_value[i] );
          if (attr_name [i])  delete[] attr_name [i];
          if (attr_value[i])  delete[] attr_value[i];
        }
        FreeVectorMemory ( attr_name ,0 );
        FreeVectorMemory ( attr_value,0 );
        attr_name  = an;
        attr_value = av;
      }

      CreateCopy ( attr_name [nAttributes],name  );
      CreateCopy ( attr_value[nAttributes],value );
      nAttributes++;

    }

    void  XMLObject::AddAttribute ( cpstr name, const int iV )  {
    char S[100];
      sprintf ( S,"%i",iV );
      AddAttribute ( name,S );
    }

    void  XMLObject::AddAttribute ( cpstr name, const bool bV )  {
      if (bV)  AddAttribute ( name,"Yes" );
         else  AddAttribute ( name,"No"  );
    }

    void  XMLObject::SetData ( cpstr Data )  {
    pstr p,d;
    int  n;
      // count ampersands
      p = pstr(Data);
      n = 0;
      while (*p)  {
        if (*p=='&')  n += 4;
        p++;
      }
      // calculate the Data length
      n += strlen(Data) + 1;  // eugene
      // allocate data space
      if (objData)  delete[] objData;
      objData = new char[n];
      // copy data, preceeding ampersands with the escape
      p = pstr(Data);
      d = objData;
      while (*p)  {
        if (*p=='&')  {
          d[0] = '&';
          d[1] = 'a';
          d[2] = 'm';
          d[3] = 'p';
          d[4] = ';';
          d += 5;
        } else  {
          *d = *p;
          d++;
        }
        p++;
      }
      *d = char(0);
    }

    void  XMLObject::AddData ( cpstr Data )  {
    pstr d1,d2;
      d1      = objData;
      objData = NULL;
      SetData ( Data );
      d2      = objData;
      objData = NULL;
      CreateConcat ( objData,d1,d2 );
    }

    void  XMLObject::SetData ( const realtype V, const int length )  {
    char N[500];
      sprintf    ( N,"%-.*g",length,V );
      CreateCopy ( objData,N );
    }

    void  XMLObject::SetData ( const int iV, const int length )  {
    char N[500];
      sprintf    ( N,"%*i",length,iV );
      CreateCopy ( objData,N );
    }

    void  XMLObject::SetData ( const bool bV )  {
      if (bV)  CreateCopy ( objData,pstr("Yes") );
         else  CreateCopy ( objData,pstr("No")  );
    }


    int XMLObject::AddMMCIFCategory ( mmcif::PCategory mmCIFCat )  {
      if (mmCIFCat->GetCategoryID()==mmcif::MMCIF_Loop)
        return AddMMCIFLoop ( mmcif::PLoop(mmCIFCat) );
      if (mmCIFCat->GetCategoryID()==mmcif::MMCIF_Struct)
        return AddMMCIFStruct ( mmcif::PStruct(mmCIFCat) );
      return -1;
    }

    pstr getCCIFTag ( pstr & ccifTag, cpstr Tag )  {
      if (Tag[0]=='_')  return CreateCopCat ( ccifTag,pstr("ccif") ,Tag );
                  else  return CreateCopCat ( ccifTag,pstr("ccif_"),Tag );
    }

    int XMLObject::AddMMCIFStruct ( mmcif::PStruct mmCIFStruct )  {
    PXMLObject XMLObject1,XMLObject2;
    pstr       SName,Tag,Field, ccifTag;
    int        nTags,i,k;

      XMLObject1 = this;

      ccifTag = NULL;

      SName = mmCIFStruct->GetCategoryName();
      if (SName)  {
        if (SName[0]!=char(1))
          XMLObject1 = new XMLObject ( getCCIFTag(ccifTag,SName) );
      }

      k = 0;
      nTags = mmCIFStruct->GetNofTags();
      for (i=0;i<nTags;i++)  {
        Tag = mmCIFStruct->GetTag ( i );
        if (Tag)  {
          XMLObject2 = new XMLObject ( getCCIFTag(ccifTag,Tag) );
          Field = mmCIFStruct->GetField ( i );
          if (Field)  {
            if (Field[0]!=char(2))  XMLObject2->SetData ( Field );
                              else  XMLObject2->SetData ( &(Field[1]) );
          }
          XMLObject1->AddObject ( XMLObject2 );
          k++;
        }
      }

      if (SName)  {
        if (SName[0]!=char(1))
          AddObject ( XMLObject1 );
      }

      if (ccifTag)  delete[] ccifTag;

      return k;

    }

    int XMLObject::AddMMCIFLoop ( mmcif::PLoop mmCIFLoop )  {
    PXMLObject XMLObject1,XMLObject2,XMLObject3;
    pstr        SName,Tag,Field,ccifTag;
    int         nTags,nRows,i,j,k;

      XMLObject1 = this;

      ccifTag = NULL;

      SName = mmCIFLoop->GetCategoryName();
      if (SName)  {
        if (SName[0]!=char(1))
          XMLObject1 = new XMLObject ( getCCIFTag(ccifTag,SName) );
      }

      k = 0;
      nTags = mmCIFLoop->GetNofTags   ();
      nRows = mmCIFLoop->GetLoopLength();
      for (i=0;i<nRows;i++)  {
        XMLObject2 = new XMLObject ( pstr("row"),
                          new XMLObject(pstr("_sernum_"),i+1) );
        for (j=0;j<nTags;j++)  {
          Tag = mmCIFLoop->GetTag ( j );
          if (Tag)  {
            XMLObject3 = new XMLObject ( getCCIFTag(ccifTag,Tag) );
            Field = mmCIFLoop->GetField ( i,j );
            if (Field)  {
              if (Field[0]!=char(2))  XMLObject3->SetData ( Field );
                                else  XMLObject3->SetData ( &(Field[1]) );
            }
            XMLObject2->AddObject ( XMLObject3 );
            k++;
          }
        }
        XMLObject1->AddObject ( XMLObject2 );
      }

      if (SName)  {
        if (SName[0]!=char(1))
          AddObject ( XMLObject1 );
      }

      if (ccifTag)  delete[] ccifTag;

      return k;

    }

    int  XMLObject::AddMMCIFData ( mmcif::PData mmCIFData )  {
    mmcif::PCategory mmCIFCat;
    int              nCats,i,k,n;
      nCats = mmCIFData->GetNumberOfCategories();
      k     = 0;
      n     = 0;
      for (i=0;(i<nCats) && (n>=0);i++)  {
        mmCIFCat = mmCIFData->GetCategory ( i );
        if (mmCIFCat)  {
          if (mmCIFCat->GetCategoryID()==mmcif::MMCIF_Loop)
            n = AddMMCIFLoop ( mmcif::PLoop(mmCIFCat) );
          else if (mmCIFCat->GetCategoryID()==mmcif::MMCIF_Struct)
            n = AddMMCIFStruct ( mmcif::PStruct(mmCIFCat) );
          else
            n = -1;
          if (n>=0)  k += n;
        }
      }
      if (n<0)  return -(k+1);
      return k;
    }

    pstr  XMLObject::GetData ( cpstr Tag, int objNo )  {
    PXMLObject XMLObject;
      XMLObject = GetObject ( Tag,objNo );
      if (XMLObject)  return XMLObject->objData;
      return NULL;
    }

    XML_RC XMLObject::GetData ( pstr & Data, cpstr Tag, int objNo )  {
    PXMLObject XMLObject;
      XMLObject = GetObject ( Tag,objNo );
      if (XMLObject)  {
        Data = XMLObject->objData;
        return XMLRC_Ok;
      } else  {
        Data = NULL;
        return XMLRC_NoTag;
      }
    }

    XML_RC XMLObject::GetData ( realtype & V, cpstr Tag, int objNo )  {
    pstr   d,p;
    XML_RC rc;
      rc = GetData ( d,Tag,objNo );
      if (d)  {
        V = strtod(d,&p);
        if ((V==0.0) && (p==d))  rc = XMLRC_RFormatError;
                           else  rc = XMLRC_Ok;
      } else if (!rc)
        rc = XMLRC_NoTag;
      return rc;
    }

    XML_RC XMLObject::GetData ( int & iV, cpstr Tag, int objNo )  {
    pstr d,p;
    XML_RC rc;
      rc = GetData ( d,Tag,objNo );
      if (d)  {
        iV = mround(strtod(d,&p));
        if ((iV==0) && (p==d))  rc = XMLRC_IFormatError;
                          else  rc = XMLRC_Ok;
      } else if (!rc)
        rc = XMLRC_NoTag;
      return rc;
    }

    XML_RC XMLObject::GetData ( bool & bV, cpstr Tag, int objNo )  {
    pstr   d;
    XML_RC rc;
      rc = GetData ( d,Tag,objNo );
      if (d)  {
        if (!strcasecmp(d,"Yes"))
          bV = true;
        else {
          bV = false;
          if (strcasecmp(d,"No"))  rc = XMLRC_OFormatError;
        }
      } else if (rc==XMLRC_Ok)
        rc = XMLRC_NoTag;
      return rc;
    }

    PXMLObject XMLObject::GetObject ( cpstr Tag, int objNo )  {
    // allow for "tag1>tag2>tag3>..."
    PXMLObject XMLObject;
    pstr       p,p1;
    int        i,j,k,l;
      XMLObject = this;
      if (Tag)  {
        p = pstr(Tag);
        do  {
          p1 = p;
          l  = 0;
          while (*p1 && (*p1!='>'))  {
            p1++;
            l++;
          }
          if (l>0)  {
            k = -1;
            j = 0;
            for (i=0;(i<XMLObject->nObjects) && (k<0);i++)
              if (XMLObject->object[i])  {
                if (!strncmp(XMLObject->object[i]->objTag,p,l))  {
                  j++;
                  if (j==objNo)  k = i;
                }
              }
            if (k<0)  {
              XMLObject = NULL;
              l = 0;
            } else  {
              XMLObject = XMLObject->object[k];
              if (*p1)  p = p1 + 1;
                  else  l = 0;
            }
          }
        } while (l>0);
      }
      return XMLObject;
    }

    PXMLObject XMLObject::GetFirstObject()  {
      if (nObjects>0)  return object[0];
      return NULL;
    }

    PXMLObject XMLObject::GetLastObject()  {
      if (nObjects>0)  return object[nObjects-1];
      return NULL;
    }

    PXMLObject XMLObject::GetObject ( int objectNo )  {
      if ((0<=objectNo) && (objectNo<nObjects))
        return object[objectNo];
      return NULL;
    }

    void  XMLObject::AddObject ( PXMLObject XMLObject, int lenInc )  {
    PPXMLObject obj1;
    int         i;

      if (!XMLObject)  return;

      if (nObjects>=nAlloc)  {
        nAlloc += lenInc;
        obj1 = new PXMLObject[nAlloc];
        for (i=0;i<nObjects;i++)
          obj1[i] = object[i];
        for (i=nObjects;i<nAlloc;i++)
          obj1[i] = NULL;
        if (object)  delete[] object;
        object = obj1;
      }

      if (object[nObjects])  delete object[nObjects];
      object[nObjects] = XMLObject;
      XMLObject->SetParent ( this );
      nObjects++;

    }


    void  XMLObject::InsertObject ( PXMLObject XMLObject, int pos,
                                    int lenInc )  {
    PPXMLObject obj1;
    int         i;

      if (!XMLObject)  return;
      if (pos>=nObjects)  {
        AddObject ( XMLObject,lenInc );
        return;
      }

      if (nObjects>=nAlloc)  {
        nAlloc += lenInc;
        obj1 = new PXMLObject[nAlloc];
        for (i=0;i<nObjects;i++)
          obj1[i] = object[i];
        for (i=nObjects;i<nAlloc;i++)
          obj1[i] = NULL;
        if (object)  delete[] object;
        object = obj1;
      }

      for (i=nObjects;i>pos;i--)
        object[i] = object[i-1];

      object[pos] = XMLObject;
      XMLObject->SetParent ( this );
      nObjects++;

    }

    XML_RC XMLObject::WriteObject ( cpstr FName, int pos, int indent )  {
    io::File f;
      f.assign ( FName,true );
      if (f.rewrite())  {
        WriteObject ( f,pos,indent );
        f.shut();
        return XMLRC_Ok;
      }
      return XMLRC_CantOpenFile;
    }

    void  XMLObject::WriteObject ( io::RFile f, int pos, int indent )  {
    int   i,pos1,lm,rm,tl;
    pstr  indstr,p,p1,q;
    bool  sngline;

      if (objTag)  {

        pos1   = pos + indent;
        indstr = new char[pos1+1];
        for (i=0;i<pos1;i++)  indstr[i] = ' ';
        indstr[pos1] = char(0);

        indstr[pos]  = char(0);
        f.Write ( indstr    );
        f.Write ( pstr("<") );
        f.Write ( objTag    );
        for (i=0;i<nAttributes;i++)  {
          f.Write ( " " );
          f.Write ( attr_name[i] );
          f.Write ( "=\"" );
          f.Write ( attr_value[i] );
          f.Write ( "\"" );
        }
        if ((!objData) && (nObjects<=0))  {
          f.WriteLine ( pstr("/>") );
          delete[] indstr;
          return;
        }
        f.Write ( pstr(">") );

        sngline = false;
        if (objData)  {
          rm = 72;               // right margin
          lm = IMin ( pos1,36 ); // left margin
          tl = strlen(objTag);
          if ((pos+tl+2+(int)strlen(objData)<rm-tl-2) &&
              (nObjects<=0))  {
            // single-line output
            sngline = true;
            f.Write ( objData );
          } else  {
            // multiple-line output with indentation
            indstr[pos] = ' ';
            indstr[lm]  = char(0);
            f.LF();
            p = objData;
            do  {
              p1 = p;
              i  = lm;
              q  = NULL;
              while ((*p1) && ((i<rm) || (!q)))  {
                if (*p1==' ')  q = p1;
                p1++;
                i++;
              }
              f.Write ( indstr );
              if (*p1)  {  // wrap data
                *q = char(0);
                f.WriteLine ( p );
                *q = ' ';
                p  = q;
                while (*p==' ')  p++;
                if (*p==char(0))  p = NULL;

              } else {  // data exchausted
                f.WriteLine ( p );
                p = NULL;
              }
            } while (p);
            indstr[lm]  = ' ';
            indstr[pos] = char(0);
          }
        } else
          f.LF();

        for (i=0;i<nObjects;i++)
          if (object[i])
            object[i]->WriteObject ( f,pos+indent,indent );

        if (!sngline)  f.Write ( indstr );
        f.Write ( pstr("</") );
        f.Write ( objTag  );
        f.WriteLine ( pstr(">") );

        delete[] indstr;

      }

    }


    XML_RC XMLObject::ReadObject ( cpstr FName )  {
    io::File f;
    char     S[500];
    int      i;
    XML_RC   rc;

      f.assign ( FName,true );
      if (f.reset(true))  {
        S[0] = char(0);
        i    = 0;
        rc   = ReadObject ( f,S,i,sizeof(S) );
        f.shut();
      } else
        rc = XMLRC_NoFile;

      if (rc!=XMLRC_Ok)  FreeMemory();

      return rc;

    }

    XML_RC XMLObject::ReadObject ( io::RFile f, pstr S,
                                   int & pos, int slen )  {
    PXMLObject xmlObject;
    pstr       S1;
    int        k,k1,k2;
    XML_RC     rc;
    bool       Done;

      FreeMemory();

      rc = XMLRC_Ok;

      k1 = -1;
      k2 = -1;
      while ((!f.FileEnd()) && (k1<0))  {
        k = strlen(S);
        while ((pos<k) && (k1<0))
          if (S[pos]=='<')  {
            if (S[pos+1]=='?') // in this version, ignore <?xxx ?>
                               //   constructions
                 pos++;
            else if (S[pos+1]!='<')
                 k1 = pos;
            else pos += 2;
          } else
            pos++;
        if (k1>=0)  {
          k2 = -1;
          while ((pos<k) && (k2<0))
            if (S[pos]=='>')  {
              if (S[pos+1]!='>')  k2 = pos;
                            else  pos += 2;
            } else
              pos++;
          if (k2<0)  rc = XMLRC_BrokenTag;
        }
        if (k1<0)  {
          f.ReadLine ( S,slen );
          pos = 0;
        }
      }

      if (k1<0)         return XMLRC_NoTag;
      if (rc!=XMLRC_Ok) return rc;

      pos++;
      if (S[k2-1]=='/')  {  // <Tag/>
        S[k2-1] = char(0);
        CreateCopy ( objTag,&(S[k1+1]) );
        return XMLRC_Ok;
      }

      S[k2] = char(0);
      CreateCopy ( objTag,&(S[k1+1]) );
      S[k2] = '>';

      S1   = new char[slen+1];
      Done = false;
      while ((!f.FileEnd()) && (!Done))  {
        k = strlen(S);
        while ((pos<k) && (!Done))  {
          k1 = pos;
          k2 = -1;
          while ((pos<k) && (k2<0))
            if (S[pos]=='<')  {
              if (S[pos+1]!='<')  k2 = pos;
                            else  pos +=2;
            } else
              pos++;
          if (k2>=0)  S[k2] = char(0);
          strcpy_des   ( S1,&(S[k1]) );
          if (S1[0])  {
            if (objData)  CreateConcat ( objData,pstr(" "),S1 );
                    else  CreateConcat ( objData,S1 );
          }
          if (k2>=0)  {
            S[k2] = '<';
            if (S[k2+1]!='/')  {
              xmlObject = new XMLObject();
              AddObject ( xmlObject );
              rc   = xmlObject->ReadObject ( f,S,pos,slen );
              Done = (rc!=XMLRC_Ok);
            } else  {
              Done = true;
              k1   = k2+2;
              k2   = -1;
              while ((pos<k) && (k2<0))
                if (S[pos]=='>')  {
                  if (S[pos+1]!='>')  k2 = pos;
                                else  pos += 2;
                  } else
                    pos++;
              if (k2<0)
                rc = XMLRC_BrokenTag;
              else  {
                S[k2] = char(0);
                if (strcmp(objTag,&(S[k1])))  rc = XMLRC_UnclosedTag;
                                        else  pos++;
              }
            }
          }
        }
        if (!Done)  {
          f.ReadLine ( S,slen );
          pos = 0;
        }
      }

      delete[] S1;

      // this keeps pairs <tag></tag> instead of replacing them for <tag/>
      // on output
      if ((!objData) && (nObjects<=0))
        CreateCopy ( objData,pstr("") );

      if (rc!=XMLRC_Ok)  FreeMemory();
      return rc;

    }


    void  XMLObject::Copy ( PXMLObject xmlObject )  {
    int i;

      FreeMemory();

      CreateCopy ( objTag ,xmlObject->objTag  );
      CreateCopy ( objData,xmlObject->objData );

      nObjects = xmlObject->nObjects;
      nAlloc   = nObjects;
      if (nObjects>0)  {
        object = new PXMLObject[nObjects];
        for (i=0;i<nObjects;i++)
          if (xmlObject->object[i])  {
            object[i] = new XMLObject();
            object[i]->Copy ( xmlObject->object[i] );
          } else
            object[i] = NULL;
      }

      nAttributes = xmlObject->nAttributes;
      nAttrAlloc  = nAttributes;
      if (nAttributes>0)  {
        GetVectorMemory ( attr_name ,nAttrAlloc,0 );
        GetVectorMemory ( attr_value,nAttrAlloc,0 );
        for (i=0;i<nAttributes;i++)  {
          attr_name [i] = NULL;
          attr_value[i] = NULL;
          CreateCopy ( attr_name [i],xmlObject->attr_name [i] );
          CreateCopy ( attr_value[i],xmlObject->attr_value[i] );
        }
      }

    }


    void  XMLObject::write ( io::RFile f )  {
    int i;
      f.CreateWrite ( objTag       );
      f.CreateWrite ( objData      );
      f.WriteInt    ( &nObjects    );
      for (i=0;i<nObjects;i++)
        StreamWrite ( f,object[i] );
      f.WriteInt    ( &nAttributes );
      for (i=0;i<nAttributes;i++)  {
        f.CreateWrite ( attr_name [i] );
        f.CreateWrite ( attr_value[i] );
      }
    }

    void  XMLObject::read ( io::RFile f )  {
    int i;

      FreeMemory();

      f.CreateRead ( objTag    );
      f.CreateRead ( objData   );

      f.ReadInt    ( &nObjects );
      nAlloc = nObjects;
      if (nObjects>0)  {
        object = new PXMLObject[nObjects];
        for (i=0;i<nObjects;i++)  {
          object[i] = NULL;
          StreamRead ( f,object[i] );
        }
      }

      f.ReadInt    ( &nAttributes );
      nAttrAlloc = nAttributes;
      if (nAttributes>0)  {
        GetVectorMemory ( attr_name ,nAttrAlloc,0 );
        GetVectorMemory ( attr_value,nAttrAlloc,0 );
        for (i=0;i<nAttributes;i++)  {
          attr_name [i] = NULL;
          attr_value[i] = NULL;
          f.CreateRead ( attr_name [i] );
          f.CreateRead ( attr_value[i] );
        }
      }

    }


    MakeStreamFunctions(XMLObject)



    PXMLObject mmCIF2XML ( mmcif::PData mmCIFData, int * rc )  {
    PXMLObject xmlObject;
    pstr       dataName;
    int        k;
      xmlObject = NULL;
      if (rc) *rc = -2;
      if (mmCIFData)  {
        dataName = mmCIFData->GetDataName();
        if (dataName)  {
          if (dataName[0])
            xmlObject = new XMLObject ( dataName );
        }
        if (!xmlObject)
          xmlObject = new XMLObject ( pstr("no_data_name") );
        k = xmlObject->AddMMCIFData ( mmCIFData );
        if (rc)  *rc = k;
      }
      return xmlObject;
    }

    PXMLObject mmCIF2XML ( cpstr XMLName, mmcif::PFile mmCIFFile,
                           int * rc )  {
    PXMLObject   xmlObject1,xmlObject2;
    mmcif::PData mmCIFData;
    int          nData,i,k,rc1;
      xmlObject1 = new XMLObject ( XMLName );
      if (rc) *rc = -1;
      if (mmCIFFile)  {
        nData = mmCIFFile->GetNofData();
        k   = 0;
        rc1 = 0;
        for (i=0;(i<nData) && (rc1>=0);i++)  {
          mmCIFData = mmCIFFile->GetCIFData ( i );
          if (mmCIFData)  {
            xmlObject2 = mmCIF2XML ( mmCIFData,&rc1 );
            if (xmlObject2)  {
              if (rc1>=0)  {
                xmlObject1->AddObject ( xmlObject2 );
                k += rc1;
              } else
                delete xmlObject2;
            }
          }
        }
        if (rc1<0)  {
          delete xmlObject1;
          if (rc)  *rc = -2;
        } else if (rc)
          *rc = k;
      }
      return xmlObject1;
    }


  }  // namespace xml

}  // namespace mmdb
