//  $Id: mmdb_mmcif_.cpp $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2015.
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
//    28.09.15   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  MMDB_MMCIF <implementation>
//       ~~~~~~~~~
//  **** Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::mmcif::Category ( mmCIF category    )
//       ~~~~~~~~~  mmdb::mmcif::Struct   ( mmCIF structure   )
//                  mmdb::mmcif::Loop     ( mmCIF loop        )
//                  mmdb::mmcif::Data     ( mmCIF data block  )
//                  mmdb::mmcif::File     ( mmCIF file        )
//
//  (C) E. Krissinel 2000-2015
//
//  =================================================================
//

#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "mmdb_mmcif_.h"

namespace mmdb  {

  namespace mmcif  {

    //  ======================  SortTags  ===============================

    void  SortTags ( psvector tag, int len, ivector index )  {
    int  i,k,l,l1,l2;
      if (len==1)  {
        index[0] = 0;
        return;
      }
      if (strcasecmp(tag[0],tag[1])<0)  {
        index[0] = 0;
        index[1] = 1;
      } else  {
        index[0] = 1;
        index[1] = 0;
      }
      for (k=2;k<len;k++)  {
        l2 = k-1;
        if (strcasecmp(tag[k],tag[index[0]])<0)       l2 = 0;
        else if (strcasecmp(tag[k],tag[index[l2]])>0) l2 = k;
        else  {
          l1 = 0;
          while (l1<l2-1)  {
            l = (l1+l2)/2;
            if (strcasecmp(tag[k],tag[index[l]])<0)  l2 = l;
                                               else  l1 = l;
          }
        }
        for (i=k;i>l2;i--)
          index[i] = index[i-1];
        index[l2] = k;
      }
    }


    //  ======================  Category  ==========================

    const int CIF_NODATA_DOT        = 0;
    const int CIF_NODATA_QUESTION   = 1;
    cpstr CIF_NODATA_DOT_FIELD      = pstr("\x02" ".");
    cpstr CIF_NODATA_QUESTION_FIELD = pstr("\x02" "?");

    Category::Category() : io::Stream()  {
      InitCategory();
    }

    Category::Category ( cpstr N ) : io::Stream()  {
      InitCategory();
      SetCategoryName ( N );
    }

    Category::Category ( io::RPStream Object ) : io::Stream(Object)  {
      InitCategory();
    }

    Category::~Category()  {
      FreeMemory();
    }

    void Category::InitCategory()  {
      name       = NULL;
      nTags      = 0;
      tag        = NULL;
      index      = NULL;
      nAllocTags = 0;
    }

    void Category::FreeMemory()  {
    int i;
      if (name)  delete[] name;
      name = NULL;
      for (i=0;i<nAllocTags;i++)
        if (tag[i])  delete[] tag[i];
      FreeVectorMemory ( tag  ,0 );
      FreeVectorMemory ( index,0 );
      nTags      = 0;
      nAllocTags = 0;
    }

    void Category::SetCategoryName ( cpstr N )  {
      if (N[0])  CreateCopy ( name,N );
      else  {
        CreateCopy ( name,pstr(" ") );
        name[0] = char(1);  // no category name
      }
    }

    void Category::ExpandTags ( int nTagsNew )  {
    int      i,nAT;
    psvector tag1;
    ivector  index1;
      if (nTagsNew>nAllocTags)  {
        nAT = nTagsNew + IMin(nAllocTags/2+1,20);
        GetVectorMemory ( tag1  ,nAT,0 );
        GetVectorMemory ( index1,nAT,0 );
        for (i=0;i<nAllocTags;i++)  {
          tag1  [i] = tag  [i];
          index1[i] = index[i];
        }
        for (i=nAllocTags;i<nAT;i++)  {
          tag1  [i] = NULL;
          index1[i] = i;
        }
        FreeVectorMemory ( tag  ,0 );
        FreeVectorMemory ( index,0 );
        tag        = tag1;
        index      = index1;
        nAllocTags = nAT;
      }
    }

    pstr Category::GetTag ( int tagNo )  {
      if ((tagNo>=0) && (tagNo<nTags))  return tag[tagNo];
      return NULL;
    }

    void Category::Sort()  {
    //  Sorts tags for easing the search
    int i,k;
      if (nAllocTags>0)  {
        k = 0;
        if (!index)
          GetVectorMemory ( index,nAllocTags,0 );
        for (i=0;i<nTags;i++)
          if (tag[i])  {
            if (k<i)  {
              tag[k] = tag[i];
              tag[i] = NULL;
            }
            k++;
          }
        nTags = k;
        SortTags ( tag,nTags,index );
      }
    }

    void Category::Optimize()  {
    int      i,k;
    psvector tag1;
      k = 0;
      for (i=0;i<nTags;i++)
        if (tag[i])  k++;
      if (k<=0)  FreeMemory();
      else if (k!=nAllocTags)  {
        GetVectorMemory  ( tag1,k,0 );
        FreeVectorMemory ( index,0 );
        k = 0;
        for (i=0;i<nTags;i++)
          if (tag[i])
            tag1[k++] = tag[i];
        FreeVectorMemory ( tag,0 );
        tag        = tag1;
        nTags      = k;
        nAllocTags = nTags;
        Sort();
      }
    }

    int  Category::GetTagNo ( cpstr ttag )  {
    //   Binary search for index of tag ttag in tag[].
    // Return:
    //    >=0 : position of the tag found
    //     <0 : the tag was not found, it could be inserted before
    //          (-RC-1)th element, where RC is the return value
    int l1,l2,l,k;

      if (!tag)    return -1;

      if (!index)  Sort();

      l = 0;
      l1 = 0;
      l2 = nTags-1;
      k  = 1;
      while (l1<l2-1)  {
        l = (l1+l2)/2;
        k = strcasecmp ( ttag,tag[index[l]] );
        if (k<0)      l2 = l;
        else if (k>0) l1 = l;
        else {
          l1 = l;
          break;
        }
      }

      if (k==0)  return index[l];    // is at RCth position
      k = strcasecmp ( ttag,tag[index[l1]] );
      if (k==0)  return index[l1];   // is at RCth position
      if (k<0)   return -1;          // would be at (-RC-1)th position
      if (l2!=l1)  {
        k = strcasecmp ( ttag,tag[index[l2]] );
        if (k==0)  return index[l2]; // is at RCth position
        if (k>0)   return -2-l2;     // would be at l2+1=(-RC-1)th position
      }

      return -2-l1;                  // would be at l1+1=(-RC-1)th position

    }

    int  Category::AddTag ( cpstr ttag )  {
    //  return -1: the tag has been added on the top of array;
    //             index is added and sorted automatically
    //        >=0: the tag is already in the array -- its position
    //             is returned
    int  i1,i;
      if (!tag)  {
        ExpandTags ( 3 );  // get space for first 3 tags
        CreateCopy ( tag[0],ttag );
        nTags  = 1;
        return -nTags;  // the tag has been added on the top of array
      }
      i1 = GetTagNo ( ttag );
      if (i1>=0)  return i1;  // non-negative returns mean that
                              // the tag is already in the array
      i1 = -i1-1;  // otherwise the tag has to be added and indexed at here
      // put new tag on the top of array and update index
      ExpandTags ( nTags+1 );
      CreateCopy ( tag[nTags],ttag );
      for (i=nTags;i>i1;i--)
        index[i] = index[i-1];
      index[i1] = nTags;
      nTags++;
      return -nTags; // the tag has been added on the top of array
    }

    void Category::PrintTags()  {
    int i;
      Sort();
      printf ( " Unsorted tags:\n" );
      for (i=0;i<nTags;i++)
        if (tag[i])
          printf ( "  %s.%s\n",name,tag[i] );
      if (index)  {
        printf ( " Sorted tags:\n" );
        for (i=0;i<nTags;i++)
          if (tag[index[i]])
            printf ( "  %s.%s\n",name,tag[index[i]] );
      }
    }

    bool Category::CheckTags ( cpstr * tagList )  {
    int i;
      i = 0;
      while (tagList[i][0])  {
        if (GetTagNo(tagList[i])<0)  return false;
        i++;
      }
      return true;
    }

    void  Category::PutCategoryName ( cpstr newName )  {
      CreateCopy ( name,newName );
    }

    void  Category::Copy ( PCategory Category )  {
    int i;
      FreeMemory();
      if (Category)  {
        CreateCopy ( name,Category->name );
        nTags      = Category->nTags;
        nAllocTags = nTags;
        if (nTags>0) {
          GetVectorMemory ( tag  ,nAllocTags,0 );
          GetVectorMemory ( index,nAllocTags,0 );
          for (i=0;i<nTags;i++)  {
            tag[i]   = NULL;
            CreateCopy ( tag[i],Category->tag[i] );
            index[i] = Category->index[i];
          }
        }
      }
    }


    void Category::write ( io::RFile f )  {
    int i;
      if (!index)  Sort();
      f.CreateWrite ( name   );
      f.WriteInt    ( &nTags );
      for (i=0;i<nTags;i++)
        f.CreateWrite ( tag[i] );
      f.WriteVector ( index,nTags,0 );
    }

    void Category::read ( io::RFile f )  {
    int i;
      FreeMemory   ();
      f.CreateRead ( name   );
      f.ReadInt    ( &nTags );
      nAllocTags = nTags;
      if (nTags>0) {
        GetVectorMemory ( tag,nTags,0 );
        for (i=0;i<nTags;i++)  {
          tag[i] = NULL;
          f.CreateRead ( tag[i] );
        }
      }
      f.CreateReadVector ( index,0 );
    }

    MakeStreamFunctions(Category)



    //  ======================  Struct  ===========================


    Struct::Struct() : Category() {
      InitStruct();
    }

    Struct::Struct ( cpstr N ) : Category(N) {
      InitStruct();
    }

    Struct::Struct ( io::RPStream Object ) : Category(Object)  {
      InitStruct();
    }

    Struct::~Struct()  {
      FreeMemory();
    }

    void Struct::FreeMemory()  {
    int i;
      for (i=0;i<nAllocTags;i++)
        if (field[i]) delete[] field[i];
      FreeVectorMemory ( field,0 );
      Category::FreeMemory();
    }

    void Struct::InitStruct()  {
      field = NULL;
    }

    void Struct::Optimize()  {
    int      i,k;
    psvector f1;
      k = 0;
      for (i=0;i<nTags;i++)
        if (!tag[i])  {
          if (field[i])  delete[] field[i];
          field[i] = NULL;
        } else if (!field[i])  {
          delete[] tag[i];
          tag[i] = NULL;
        } else
          k++;
      if (k<=0)  FreeMemory();
      else if (k!=nAllocTags)  {
        f1 = new pstr[k];
        k  = 0;
        for (i=0;i<nTags;i++)
          if (tag[i])
            f1[k++] = field[i];
        FreeVectorMemory ( field,0 );
        field = f1;
        Category::Optimize();
      }
    }

    void Struct::AddField ( cpstr F, cpstr T, bool Concatenate )  {
    psvector field1;
    int      i,nAT;
    pstr     nf;

      nAT = nAllocTags;
      i   = AddTag ( T );

      if (i<0) {
        // The tag was not in the list, but has been added on the top
        // of list. Now expand the field list and put new field on
        // the top of it.
        if (nAllocTags>nAT)  {
          GetVectorMemory ( field1,nAllocTags,0 );
          for (i=0;i<nTags-1;i++)
            field1[i] = field[i];
          for (i=nTags-1;i<nAllocTags;i++)
            field1[i] = NULL;
          FreeVectorMemory ( field,0 );
          field = field1;
        }
        i        = nTags-1;
        field[i] = NULL;
      }

      if (!F)  {
        if ((!Concatenate) || (!field[i]))  {
          CreateCopy ( field[i],pstr(" ?") );
          field[i][0] = char(2);
        }
      } else if ((!Concatenate) || (!field[i]))
        CreateCopy ( field[i],F );
      else  {
        nf = new char[strlen(field[i])+strlen(F)+1];
        strcpy ( nf,field[i] );
        strcat ( nf,F        );
        delete[] field[i];
        field[i] = nf;
      }

    }

    pstr Struct::GetField ( int tagNo )  {
      if ((tagNo>=0) && (tagNo<nTags))  return field[tagNo];
      return NULL;
    }

    int  Struct::GetString ( pstr & S, cpstr TName,
                                   bool Remove )  {
    int k = GetTagNo ( TName );
      if (S)  delete[] S;
      S = NULL;
      if (!field)     return CIFRC_NoField;
      if (k<0)        return CIFRC_NoTag;
      if (!field[k])  return CIFRC_NoField;
      if (field[k][0]==char(2))  {
        if (Remove)  {
          delete[] field[k];
          field[k] = NULL;
        }
      } else if (Remove)  {
        S = field[k];
        field[k] = NULL;
      } else
        CreateCopy ( S,field[k] );
      return 0;
    }

    pstr Struct::GetString ( cpstr TName, int & RC )  {
    int k = GetTagNo ( TName );
      if (k<0)  {
        RC = CIFRC_NoTag;
        return NULL;
      }
      if (!field)  {
        RC = CIFRC_NoField;
        return NULL;
      }
      if (!field[k])  {
        RC = CIFRC_NoField;
        return NULL;
      }
      RC = 0;
      if (field[k][0]==char(2))  return NULL;
      return field[k];
    }

    int  Struct::DeleteField ( cpstr TName )  {
    int k = GetTagNo ( TName );
      if ((k>=0) && (field)) {
        if (field[k])  delete[] field[k];
        field[k] = NULL;
      }
      return k;
    }

    int  Struct::GetReal ( realtype & R, cpstr TName,
                                 bool Remove )  {
    pstr endptr;
    int  RC;
    int  k = GetTagNo ( TName );
      R = 0.0;
      if (!field)                return CIFRC_NoField;
      if (k<0)                   return CIFRC_NoTag;
      if (!field[k])             return CIFRC_NoField;
      if (field[k][0]==char(2))  return CIFRC_NoData;
      R = strtod ( field[k],&endptr );
      if (endptr==field[k])  RC = CIFRC_WrongFormat;
      else  {
        RC = 0;
        if (Remove)  {
          delete[] field[k];
          field[k] = NULL;
        }
      }
      return RC;
    }

    int  Struct::GetInteger ( int & I, cpstr TName,
                                    bool Remove )  {
    pstr endptr;
    int  RC;
    int  k = GetTagNo ( TName );
      I = 0;
      if (!field)                return CIFRC_NoField;
      if (k<0)                   return CIFRC_NoTag;
      if (!field[k])             return CIFRC_NoField;
      if (field[k][0]==char(2))  {
        if (field[k][1]=='.')  I = MinInt4;
        return CIFRC_NoData;
      }
      I = mround ( strtod(field[k],&endptr) );
      if (endptr==field[k])  RC = CIFRC_WrongFormat;
      else  {
        RC = 0;
        if (Remove)  {
          delete[] field[k];
          field[k] = NULL;
        }
      }
      return RC;
    }


    void Struct::PutString ( cpstr S, cpstr T,
                                   bool NonBlankOnly )  {
    pstr p;
      if (!S)  PutNoData ( CIF_NODATA_QUESTION,T );
      else  {
        p = pstr(S);
        if (NonBlankOnly)
          while (*p==' ')  p++;
        if (!(*p))  PutNoData ( CIF_NODATA_DOT,T );
              else  AddField  ( S,T,false );
      }
    }

    void Struct::PutDate ( cpstr T )  {
    time_t t;
    tm *   tstruct;
     char   S[100];
      t       = time ( NULL );
      tstruct = localtime(&t);
      if (tstruct)
            sprintf ( S,"%4i-%02i-%02i",
                tstruct->tm_year+1900,tstruct->tm_mon+1,tstruct->tm_mday );
      else  strcpy  ( S,"YYYY-MM-DD" );
      AddField ( S,T,false );
    }


    void Struct::PutNoData ( int NoDataType, cpstr T )  {
    char S[10];
      S[0] = char(2);
      if (NoDataType==CIF_NODATA_DOT)  S[1] = '.';
                                 else  S[1] = '?';
      S[2] = char(0);
      AddField ( S,T,false );
    }


    void Struct::PutReal ( realtype R, cpstr T, int prec )  {
    char rS[100];
      sprintf  ( rS,"%.*g",prec,R );
      AddField ( rS,T,false );
    }

    void Struct::PutReal ( realtype R, cpstr T, cpstr format )  {
    char rS[100];
      sprintf  ( rS,format,R );
      AddField ( DelSpaces(rS,' '),T,false );
    }

    void Struct::PutInteger ( int I, cpstr T )  {
    char iS[100];
      if (I>MinInt4)  {
        sprintf  ( iS,"%i",I );
        AddField ( iS,T,false );
      } else
        PutNoData ( CIF_NODATA_DOT,T );
    }



    #define  NODATA_Q  pstr("?")
    #define  NODATA_P  pstr(".")

    bool Struct::WriteMMCIFStruct ( cpstr FName,
                                    io::GZ_MODE gzipMode )  {
    io::File f;
      f.assign ( FName,true,false,gzipMode );
      if (f.rewrite())  {
        WriteMMCIF ( f );
        f.shut();
        return true;
      } else
        return false;
    }

    #define _max_output_line_width 256

    void Struct::WriteMMCIF ( io::RFile f )  {
    int   i,j,k,l,m,n;
    pstr  F;

      // calculate maximal length of tags
      l = 0;
      for (i=0;i<nTags;i++)
        l = IMax(l,strlen(tag[i]));
      l += 1;  // add one space separator

      // calculate maximal space left for data
      m  = _max_output_line_width - l;
      // subtract category name width
      if (name[0]!=char(1))  m -= strlen(name);

      // start outout
      f.LF();
      for (i=0;i<nTags;i++)  {  // for each tag

        // print category name, if not hidden, and dot
        if (name[0]!=char(1))  {
          f.Write ( name      );
          f.Write ( pstr(".") );
        }

        // print tag, checking for duplicate tag flag
        F = FirstOccurence ( tag[i],'\1' );
        if (F)  {
          *F = char(0);
          f.Write ( tag[i] );
          *F = '\1';
        } else
          f.Write ( tag[i] );

        // print field
        if (field[i])  {  // field is defined
          F = field[i];
          if (FirstOccurence(F,'\n') || strstr(F,"\" "))  {
            f.Write ( pstr("\n;")   );
            f.Write ( F             );
            f.Write ( pstr("\n;\n") );
          } else {
            n = strlen(F);
            if (n>m)  // wrap around if field is too long
              f.Write ( pstr("\n ") );
            else {
              k = l-strlen(tag[i]);
              for (j=0;j<k;j++)
                f.Write ( pstr(" ") );
            }
            if ((((F[0]=='.') || (F[0]=='?')) && (!F[1])) ||
                FirstOccurence(F,' '))  {
              f.Write ( pstr("\"")   );
              f.Write ( field[i]     );
              f.Write ( pstr("\"\n") );
            } else if (field[i][0]==char(2))  {
              f.WriteLine ( &(field[i][1]) );
            } else if (!field[i][0])  {
              f.WriteLine ( NODATA_P );
            } else
              f.WriteLine ( field[i] );
          }

        } else  {  // field if not defined, put question mark

          k = l-strlen(tag[i]);
          for (j=0;j<k;j++)
            f.Write ( pstr(" ") );
          f.WriteLine ( NODATA_Q );

        }

      }

    }

    void Struct::Copy ( PCategory Struct )  {
    int i;
      Category::Copy ( Struct );
      if (nTags>0)  {
        GetVectorMemory ( field,nTags,0 );
        for (i=0;i<nTags;i++)  {
          field[i] = NULL;
          CreateCopy ( field[i],PStruct(Struct)->field[i] );
        }
      }
    }

    void Struct::write ( io::RFile f )  {
    int i;
      Category::write ( f );
      for (i=0;i<nTags;i++)
        f.CreateWrite ( field[i] );
    }

    void Struct::read ( io::RFile f )  {
    int i;
      Category::read ( f );
      if (nTags>0)  {
        GetVectorMemory ( field,nTags,0 );
        for (i=0;i<nTags;i++)  {
          field[i] = NULL;
          f.CreateRead ( field[i] );
        }
      }
    }


    MakeStreamFunctions(Struct)



    //  ======================  Loop  ==============================


    Loop::Loop() : Category()  {
      InitLoop();
    }

    Loop::Loop ( cpstr N ) : Category(N)  {
      InitLoop();
    }

    Loop::Loop ( io::RPStream Object ) : Category(Object)  {
      InitLoop();
    }

    Loop::~Loop()  {
      FreeMemory();
    }

    void Loop::InitLoop()  {
      nRows      = 0;
      field      = NULL;
      iColumn    = 0;
      nAllocRows = 0;
    }

    void Loop::FreeMemory()  {
      DeleteFields();
      Category::FreeMemory();
    }

    void Loop::Optimize()  {
    int      i,j,nT,nR,k,m;
    bool  empty;
    psmatrix f1;

      if (!field)  {
        Category::Optimize();  // optimize tags
        return;
      }

      // first check for empty columns
      nT = 0;
      for (i=0;i<nTags;i++)
        if (!tag[i])  {
          for (j=0;j<nRows;j++)  // delete ith column of field
            if (field[j])  {
              if (field[j][i])
                delete[] field[j][i];
              field[j][i] = NULL;
            }
        } else  {
          empty = true;
          j = 0;
          while ((j<nRows) && empty)  {  // check if ith column is empty
            if (field[j])
              empty = !field[j][i];
            j++;
          }
          if (empty)  {    // if ith column is empty, delete its tag
            delete[] tag[i];
            tag[i] = NULL;
          } else           // otherwise count ith tag
            nT++;
        }

      // now check for empty rows
      nR = 0;
      for (j=0;j<nRows;j++)
        if (field[j])  {
          i = 0;
          while ((i<nTags) && (!field[j][i])) i++;
          if (i>=nTags)  {
            delete[] field[j];  // delete empty row
            field[j] = NULL;
          } else
            nR++;             // count non-empty row
        }
      if ((nT<=0) || (nR<=0))
        FreeMemory();  // the loop is completely empty
      else if ((nT!=nTags) || (nR!=nAllocRows))  {
        f1 = new psvector[nR];
        m  = 0;
        for (j=0;j<nRows;j++)
          if (field[j])  {
            f1[m] = new pstr[nT];
            k = 0;
            for (i=0;i<nTags;i++)
              if (tag[i])
                f1[m][k++] = field[j][i];
            m++;
            delete[] field[j];
          }
        if (field)  delete[] field;
        field = f1;
        nRows = nR;
        nAllocRows = nRows;
        Category::Optimize();  // optimize tags
      }

    }

    void Loop::DeleteFields()  {
    int i,j;
      if (field)  {
        for (i=0;i<nAllocRows;i++)
          if (field[i])  {
            for (j=0;j<nTags;j++)
              if (field[i][j])  delete[] field[i][j];
            delete[] field[i];
          }
        delete[] field;
        field      = NULL;
        nRows      = 0;
        nAllocRows = 0;
      }
    }

    void Loop::AddLoopTag ( cpstr T, bool Remove )  {
    psmatrix  f1;
    int       i,j,nT1;
      if (Remove)  {
        DeleteFields();
        AddTag ( T );
      } else  {
        f1    = field;
        field = NULL;
        i     = AddTag ( T );
        if ((f1) && (i<0))  {
          // The tag was added on the top of tag array. Create
          // and fill new fields.
          field = new psvector[nAllocRows];
          nT1   = nTags-1;
          for (i=0;i<nAllocRows;i++)
            if (f1[i])  {
              field[i] = new pstr[nTags];
              for (j=0;j<nT1;j++)
                field[i][j] = f1[i][j];
              field[i][nT1] = NULL;
              f1[i] = NULL;
            } else
              field[i] = NULL;
          delete[] f1;
        } else
          // The tag was already in the category. Just restore fields.
          field = f1;
      }
    }


    void Loop::ExpandRows ( int nRowsNew )  {
    int      nAR,i;
    psmatrix field1;
      if (nRowsNew>nAllocRows)  {
        nAR    = nRowsNew + IMin(nAllocRows/2+10,2000);
        field1 = new psvector[nAR];
        for (i=0;i<nAllocRows;i++)
          field1[i] = field[i];
        for (i=nAllocRows;i<nAR;i++)
          field1[i] = NULL;
        if (field)  delete[] field;
        field      = field1;
        nAllocRows = nAR;
      }
    }

    void Loop::AddString ( cpstr S, bool NonBlankOnly )  {
    int  i;
    pstr p;
      if (!S)  AddNoData ( CIF_NODATA_QUESTION );
      else  {
        p = pstr(S);
        if (NonBlankOnly)
          while (*p==' ')  p++;
        if (!(*p))  AddNoData ( CIF_NODATA_DOT );
        else  {
          if (iColumn==0)  {  // start a new row
            ExpandRows ( nRows+1 );
            field[nRows] = new pstr[nTags];
            for (i=0;i<nTags;i++)
              field[nRows][i] = NULL;
            nRows++;
          }
          CreateCopy ( field[nRows-1][iColumn],S );
          iColumn++;
          if (iColumn>=nTags) iColumn = 0;
        }
      }
    }

    void Loop::AddNoData ( int NoDataType )  {
    char S[10];
      S[0] = char(2);
      if (NoDataType==CIF_NODATA_DOT)  S[1] = '.';
                                 else  S[1] = '?';
      S[2] = char(0);
      AddString ( S );
    }

    void Loop::AddReal ( realtype R, int prec )  {
    char rS[100];
      sprintf ( rS,"%.*g",prec,R );
      AddString ( rS );
    }

    void Loop::AddReal ( realtype R, cpstr format )  {
    char rS[100];
      sprintf ( rS,format,R );
      AddString ( DelSpaces(rS,' ') );
    }

    void Loop::AddInteger ( int I )  {
    char iS[100];
      if (I>MinInt4)  {
        sprintf ( iS,"%i",I );
        AddString ( iS );
      } else
        AddNoData ( CIF_NODATA_DOT );
    }


    pstr Loop::GetField ( int rowNo, int tagNo )  {
      if ((tagNo>=0) && (tagNo<nTags) &&
          (rowNo>=0) && (rowNo<nRows))  {
        if (field[rowNo])  return field[rowNo][tagNo];
      }
      return NULL;
    }

    int  Loop::GetString ( pstr & S, cpstr TName, int nrow,
                                 bool Remove)  {
    int k = GetTagNo ( TName );
      if (S)  delete[] S;
      S = NULL;
      if (k<0)                       return CIFRC_NoTag;
      if ((nrow<0) || (nrow>=nRows)) return CIFRC_WrongIndex;
      if (!field[nrow])              return CIFRC_NoField;
      if (!field[nrow][k])           return CIFRC_NoField;
      if (field[nrow][k][0]==char(2))  {
        if (Remove)  {
          delete[] field[nrow][k];
          field[nrow][k] = NULL;
        }
      } else if (Remove)  {
        S = field[nrow][k];
        field[nrow][k] = NULL;
      } else
        CreateCopy ( S,field[nrow][k] );
      return 0;
    }

    pstr Loop::GetString ( cpstr TName, int nrow, int & RC )  {
    int k = GetTagNo ( TName );
      if (k<0)  {
        RC = CIFRC_NoTag;
        return NULL;
      }
      if ((nrow<0) || (nrow>=nRows))  {
        RC = CIFRC_WrongIndex;
        return NULL;
      }
      if (!field[nrow])  {
        RC = CIFRC_NoField;
        return NULL;
      }
      if (!field[nrow][k])  {
        RC = CIFRC_NoField;
        return NULL;
      }
      RC = 0;
      // char(2) means the field was either '.' or '?'
      if (field[nrow][k][0]==char(2))  return NULL;
      return field[nrow][k];
    }

    //  CopyString() does nothing if RC is not 0
    void Loop::CopyString   ( pstr  buf,   int maxlength,
                                    cpstr TName, int nrow, int & RC )  {
    pstr p;
    int  k;

      if (RC)  return;

      k = GetTagNo ( TName );
      if (k<0)  {
        RC = CIFRC_NoTag;
        buf[0] = char(0);
        return;
      }
      if ((nrow<0) || (nrow>=nRows))  {
        RC = CIFRC_WrongIndex;
        buf[0] = char(0);
        return;
      }
      if (!field[nrow])  {
        RC = CIFRC_NoField;
        buf[0] = char(0);
        return;
      }
      p = field[nrow][k];
      if (!p)  {
        RC = CIFRC_NoField;
        buf[0] = char(0);
        return;
      }

      // char(2) means the field was either '.' or '?'
      if (p[0]==char(2))  {
        buf[0] = p[0];
        buf[1] = char(0);
      } else
        strncpy ( buf,p,IMin(maxlength,strlen(p)+1) );

    }



    int Loop::DeleteField ( cpstr TName, int nrow )  {
    int k = GetTagNo ( TName );
      if (k<0)  return CIFRC_NoTag;
      if ((nrow<0) || (nrow>=nRows))
                return CIFRC_WrongIndex;
      if (field[nrow])  {
        if (field[nrow][k])  delete[] field[nrow][k];
        field[nrow][k] = NULL;
      }
      return k;
    }

    int Loop::DeleteRow ( int nrow )  {
    int i;
      if ((nrow<0) || (nrow>=nRows))
                return CIFRC_WrongIndex;
      if (field[nrow])  {
        for (i=0;i<nTags;i++)
          if (field[nrow][i])  {
            delete[] field[nrow][i];
            field[nrow][i] = NULL;
          }
        delete[] field[nrow];
        field[nrow] = NULL;
      }
      return 0;
    }

    int  Loop::GetReal ( realtype & R, cpstr TName, int nrow,
                               bool Remove )  {
    pstr endptr;
    int  k = GetTagNo ( TName );
      if (k<0)  return CIFRC_NoTag;
      if ((nrow<0) || (nrow>=nRows))
                return CIFRC_WrongIndex;
      R = 0.0;
      if (!field[nrow])                return CIFRC_NoField;
      if (!field[nrow][k])             return CIFRC_NoField;
      if (field[nrow][k][0]==char(2))  return CIFRC_NoField;
      R = strtod ( field[nrow][k],&endptr );
      if (endptr==field[nrow][k])      return CIFRC_WrongFormat;
      if (Remove)  {
        delete[] field[nrow][k];
        field[nrow][k] = NULL;
      }
      return 0;
    }

    void Loop::CopyReal ( realtype & R, cpstr TName, int nrow,
                                int & RC )  {
    pstr endptr;
    int  k;

      if (RC)  return;

    //  R = 0.0;
      k = GetTagNo ( TName );

      if (k<0)                              RC = CIFRC_NoTag;
      else if ((nrow<0) || (nrow>=nRows))   RC = CIFRC_WrongIndex;
      else if (!field[nrow])                RC = CIFRC_NoField;
      else if (!field[nrow][k])             RC = CIFRC_NoField;
      else if (field[nrow][k][0]==char(2))  RC = CIFRC_NoField;
      else  {
        R = strtod ( field[nrow][k],&endptr );
        if (endptr==field[nrow][k])  RC = CIFRC_WrongFormat;
      }

    }

    void Loop::CopyInteger ( int & I, cpstr TName, int nrow,
                                   int & RC )  {
    pstr endptr;
    int  k;

      if (RC)  return;

      I = 0;
      k = GetTagNo ( TName );

      if (k<0)                              RC = CIFRC_NoTag;
      else if ((nrow<0) || (nrow>=nRows))   RC = CIFRC_WrongIndex;
      else if (!field[nrow])                RC = CIFRC_NoField;
      else if (!field[nrow][k])             RC = CIFRC_NoField;
      else if (field[nrow][k][0]==char(2))  RC = CIFRC_NoField;
      else  {
        I = mround ( strtod ( field[nrow][k],&endptr ) );
        if (endptr==field[nrow][k])  RC = CIFRC_WrongFormat;
      }

    }

    int  Loop::GetInteger ( int & I, cpstr TName, int nrow,
                                  bool Remove )  {
    pstr endptr;
    int  k = GetTagNo ( TName );
      if (k<0)  return CIFRC_NoTag;
      if ((nrow<0) || (nrow>=nRows))
                return CIFRC_WrongIndex;
      I = 0;
      if (!field[nrow])                return CIFRC_NoField;
      if (!field[nrow][k])             return CIFRC_NoField;
      if (field[nrow][k][0]==char(2))  {
        if (field[nrow][k][1]=='.')  I = MinInt4;
        return CIFRC_NoField;
      }
      I = mround ( strtod(field[nrow][k],&endptr) );
      if (endptr==field[nrow][k])      return CIFRC_WrongFormat;
      if (Remove)  {
        delete[] field[nrow][k];
        field[nrow][k] = NULL;
      }
      return 0;
    }


    int  Loop::GetSVector ( psvector & S, cpstr TName,
                                  int i1, int i2, bool Remove )  {
    int j,k,r1,r2;
      r1 = IMin(i1,i2);
      r2 = IMin(IMax(i1,i2),nRows-1);
      if ((r1<0) || (r1>=nRows) || (r2<0)) return CIFRC_WrongIndex;
      k = GetTagNo ( TName );
      if (k<0)  return CIFRC_NoTag;
      if (!S)
        GetVectorMemory ( S,r2-r1+1,r1 );
      if (Remove)  {
        for (j=r1;j<=r2;j++)
          if (field[j])  {
            S[j] = field[j][k];
            field[j][k] = NULL;
            if (S[j])  {
              if (S[j][0]==char(2))  {
                delete[] S[j];
                S[j] = NULL;
              }
            }
          } else
            S[j] = NULL;
      } else  {
        for (j=r1;j<=r2;j++)  {
          S[j] = NULL;
          if (field[j])  {
            if (field[j][k])  {
              if (field[j][k][0]!=char(2))
                CreateCopy ( S[j],field[j][k] );
            }
          }
        }
      }
      return 0;
    }

    int  Loop::GetRVector ( rvector & R, cpstr TName,
                                  int i1, int i2, bool Remove )  {
    int  j,k,r1,r2,RC;
    pstr endptr;
      r1 = IMin(i1,i2);
      r2 = IMin(IMax(i1,i2),nRows-1);
      if ((r1<0) || (r1>=nRows) || (r2<0)) return CIFRC_WrongIndex;
      k = GetTagNo ( TName );
      if (k<0)  return CIFRC_NoTag;
      if (!R)
        GetVectorMemory ( R,r2-r1+1,r1 );
      RC = 0;
      for (j=r1;j<=r2;j++)  {
        R[j] = 0.0;
        if (field[j])  {
          if (field[j][k])  {
            R[j] = strtod ( field[j][k],&endptr );
            if (endptr==field[j][k])  RC = CIFRC_WrongFormat;
            if (Remove)  {
              delete[] field[j][k];
              field[j][k] = NULL;
            }
          }
        }
      }
      return RC;
    }

    int  Loop::GetIVector ( ivector & I, cpstr TName,
                                  int i1, int i2, bool Remove )  {
    int  j,k,r1,r2,RC;
    pstr endptr;
      r1 = IMin(i1,i2);
      r2 = IMin(IMax(i1,i2),nRows-1);
      if ((r1<0) || (r1>=nRows) || (r2<0)) return CIFRC_WrongIndex;
      k = GetTagNo ( TName );
      if (k<0)    return CIFRC_NoTag;
      if (!I)
        GetVectorMemory ( I,r2-r1+1,r1 );
      RC = 0;
      for (j=r1;j<=r2;j++)  {
        I[j] = 0;
        if (field[j])  {
          if (field[j][k])  {
            I[j] = mround ( strtod(field[j][k],&endptr) );
            if (endptr==field[j][k]) RC = CIFRC_WrongFormat;
            if (Remove)  {
              delete[] field[j][k];
              field[j][k] = NULL;
            }
          }
        }
      }
      return RC;
    }


    void Loop::PutString ( cpstr S, cpstr T, int nrow )  {
    psmatrix field1;
    int      nT,nR,iT,i,j;
      nT = nTags;
      nR = nRows;
      iT = AddTag ( T );
      if (iT<0)  iT = nTags-1;
      if (nTags>nT)  {
        // a new tag has been added; all field must be reallocated.
        nRows      = IMax(nR,nrow+1);  // nrow is indexed like 0,1,...
        nAllocRows = IMax(nR,nrow+IMin(nR/2+1,2000));
        field1     = new psvector[nAllocRows];
        for (i=0;i<nR;i++)
          if (field[i])  {
            field1[i] = new pstr[nTags];
            for (j=0;j<nT;j++)
              field1[i][j] = field[i][j];
            for (j=nT;j<nTags;j++)
              field1[i][j] = NULL;
            delete[] field[i];
          } else
            field1[i] = NULL;
        for (i=nR;i<nRows;i++)
          field1[i] = NULL;
        if (field)  delete[] field;
        field = field1;
      } else if (nrow>=nR)  {
        // only new rows are to be added
        ExpandRows ( nrow+1 );
        nRows++;
      }
      if (!field[nrow])  {
        field[nrow] = new pstr[nTags];
        for (j=0;j<nTags;j++)
          field[nrow][j] = NULL;
      }
      CreateCopy ( field[nrow][iT],S );
      iColumn = iT+1;
      if (iColumn>=nTags) iColumn = 0;
    }


    void Loop::PutNoData ( int NoDataType, cpstr T, int nrow )  {
    char S[10];
      S[0] = char(2);
      if (NoDataType==CIF_NODATA_DOT)  S[1] = '.';
                                 else  S[1] = '?';
      S[2] = char(0);
      PutString ( S,T,nrow );
    }


    void Loop::PutReal ( realtype R, cpstr T, int nrow, int prec )  {
    char rS[100];
      sprintf ( rS,"%.*g",prec,R );
      PutString ( rS,T,nrow );
    }

    void Loop::PutReal ( realtype R, cpstr T, int nrow,
                               cpstr format )  {
    char rS[100];
      sprintf ( rS,format,R );
      PutString ( DelSpaces(rS,' '),T,nrow );
    }

    void Loop::PutInteger ( int I, cpstr T, int nrow )  {
    char iS[100];
      if (I>MinInt4)  {
        sprintf ( iS,"%i",I );
        PutString ( iS,T,nrow );
      } else
        PutNoData ( CIF_NODATA_DOT,T,nrow );
    }

    void Loop::PutSVector ( psvector S, cpstr T, int i1, int i2 )  {
    int i,j,k;
      PutString ( S[i2],T,i2 );
      if (iColumn==0)  k = nTags-1;
                 else  k = iColumn-1;
      for (i=i2-1;i>=i1;i--)  {
        if (!field[i])  {
          field[i] = new pstr[nTags];
          for (j=0;j<nTags;j++)
            field[i][j] = NULL;
        }
        CreateCopy ( field[i][k],S[i] );
      }
    }

    void Loop::PutRVector ( rvector R, cpstr T,
                                  int i1, int i2, int prec )  {
    int  i,j,k;
    char rS[100];
      PutReal ( R[i2],T,i2,prec );
      if (iColumn==0)  k = nTags-1;
                 else  k = iColumn-1;
      for (i=i2-1;i>=i1;i--)  {
        if (!field[i])  {
          field[i] = new pstr[nTags];
          for (j=0;j<nTags;j++)
            field[i][j] = NULL;
        }
        sprintf ( rS,"%.*g",prec,R[i] );
        CreateCopy ( field[i][k],rS );
      }
    }

    void Loop::PutIVector ( ivector I, cpstr T,
                                  int i1, int i2 )  {
    int  l,j,k;
    char iS[100];
      PutInteger ( I[i2],T,i2 );
      if (iColumn==0)  k = nTags-1;
                 else  k = iColumn-1;
      for (l=i2-1;l>=i1;l--)  {
        if (!field[l])  {
          field[l] = new pstr[nTags];
          for (j=0;j<nTags;j++)
            field[l][j] = NULL;
        }
        sprintf ( iS,"%i",I[l] );
        CreateCopy ( field[l][k],iS );
      }
    }

    bool Loop::WriteMMCIFLoop ( cpstr FName, io::GZ_MODE gzipMode )  {
    io::File f;
      f.assign ( FName,true,false,gzipMode );
      if (f.rewrite())  {
        WriteMMCIF ( f );
        f.shut();
        return true;
      } else
        return false;
    }

    void Loop::WriteMMCIF ( io::RFile f )  {
    int     i,j,k,m,n;
    ivector l;
    pstr    F;

      // write loop keyword
      f.Write ( pstr("\nloop_\n") );

      GetVectorMemory ( l,nTags,0 );
      k = 0;
      for (i=0;i<nTags;i++)  {
        if (name[0]!=char(1))  {
          f.Write ( name      );
          f.Write ( pstr(".") );
        }
        F = FirstOccurence ( tag[i],'\1' );
        if (F)  {
          *F = char(0);
          f.WriteLine ( tag[i] );
          *F = '\1';
        } else
          f.WriteLine ( tag[i] );
        l[i] = 0;
        for (j=0;j<nRows;j++)
          if (field[j])  {
            if (field[j][i])  {
              F = field[j][i];
              if (FirstOccurence(F,'\n') || strstr(F,"\" "))
                                       l[i] = 10001;
              else if (F[0]==char(2))  l[i] = IMax(l[i],1);
              else if (((F[0]=='.') || (F[0]=='?')) &&
                       (!F[1]))        l[i] = IMax(l[i],3);
              else  {
                if (FirstOccurence(F,' ') || FirstOccurence(F,'"') || FirstOccurence(F,'\''))
                      m = 2;
                else  m = 0;
                l[i] = IMax(l[i],strlen(F)+m);
              }
            }
          }
        l[i] = IMax(l[i],1);
        k += l[i]+1;
        if (k>_max_output_line_width)  {
          l[i] = -l[i];
          k = 0;
        }
      }
      for (i=0;i<nRows;i++)  {
        m = 0;  // counts symbols in the string
        k = 0;  // rest of left-aligned fields to fill with spaces
        for (j=0;j<nTags;j++)  {
          n = k;
          k = l[j];   // length of the field
          if (k<0)  k = -k;
          m += k+1;
          if (m>_max_output_line_width)  {
            f.LF();
            m = k+1;
          } else
            while (n>0) {
              f.Write ( pstr(" ") );
              n--;
            }
          if (field[i])  {
            if (field[i][j])  {
              F = field[i][j];
              if (k>10000)  {
                if (F[0]==char(2))  {
                  f.Write     ( pstr(" ") );
                  f.WriteLine ( &(F[1])   );
                } else if (!F[0])  {
                  f.Write     ( pstr(" ") );
                  f.WriteLine ( NODATA_P  );
                } else  {
                  f.Write     ( pstr(";") );
                  f.WriteLine ( F         );
                  f.WriteLine ( pstr(";") );
                }
                m = 0;
                k = 0;
              } else if ((((F[0]=='.') ||
                          (F[0]=='?')) && (!F[1])) ||
                         FirstOccurence(F,' ') ||
                         FirstOccurence(F,'"') ||
                         FirstOccurence(F,'\''))  {
                f.Write ( pstr(" \"") );
                f.Write ( F           );
                f.Write ( pstr("\"")  );
                k -= strlen(F)+2;
              } else if (F[0]==char(2))  {
                f.Write ( pstr(" ") );
                f.Write ( &(F[1])   );
                k--;
              } else if (!F[0])  {
                f.Write ( pstr(" ") );
                f.Write ( NODATA_P  );
                k--;
              } else  {
                f.Write ( pstr(" ") );
                f.Write ( F         );
                k -= strlen(F);
              }
            } else  {
              f.Write ( pstr(" ") );
              f.Write ( NODATA_Q  );
              k--;
            }
          } else  {
            f.Write ( pstr(" ") );
            f.Write ( NODATA_Q  );
            k--;
          }
        }
        if (m) f.LF();
      }
      f.WriteLine ( pstr("#") );

    }


    void Loop::Copy ( PCategory Loop )  {
    int i,j;
      Category::Copy ( Loop );
      nRows      = PLoop(Loop)->nRows;
      nAllocRows = nRows;
      if ((nTags>0) && (nRows>0))  {
        field = new psvector[nRows];
        for (i=0;i<nRows;i++)  {
          if (PLoop(Loop)->field[i])  {
            field[i] = new pstr[nTags];
            for (j=0;j<nTags;j++)  {
              field[i][j] = NULL;
              CreateCopy ( field[i][j],PLoop(Loop)->field[i][j] );
            }
          } else
            field[i] = NULL;
        }
      }
      iColumn = PLoop(Loop)->iColumn;
    }


    void Loop::write ( io::RFile f )  {
    int i,j;
      Category::write ( f );
      f.WriteInt ( &nRows );
      if ((nTags>0) && (nRows>0))  {
        for (i=0;i<nRows;i++)  {
          if (field[i])  {
            j = 1;
            f.WriteInt ( &j );
            for (j=0;j<nTags;j++)
              f.CreateWrite ( field[i][j] );
          } else  {
            j = 0;
            f.WriteInt ( &j );
          }
        }
      }
      f.WriteInt ( &iColumn );
    }

    void Loop::read ( io::RFile f )  {
    int i,j;
      Category::read ( f );
      f.ReadInt ( &nRows );
      nAllocRows = nRows;
      if ((nTags>0) && (nRows>0))  {
        field = new psvector[nRows];
        for (i=0;i<nRows;i++)  {
          f.ReadInt ( &j );
          if (j)  {
            field[i] = new pstr[nTags];
            for (j=0;j<nTags;j++)  {
              field[i][j] = NULL;
              f.CreateRead ( field[i][j] );
            }
          } else
            field[i] = NULL;
        }
      }
      f.ReadInt ( &iColumn );
    }


    MakeStreamFunctions(Loop)



    //  ======================  Data  =============================


    Data::Data() : io::Stream()  {
      InitData();
    }

    Data::Data ( cpstr N ) : io::Stream()  {
      InitData();
      CreateCopy ( name,N );
    }

    Data::Data ( io::RPStream Object ) : io::Stream(Object)  {
      InitData();
    }

    Data::~Data() {
      FreeMemory(0);
    }

    void Data::InitData()  {
      name         = NULL;
      nCategories  = 0;
      Category     = NULL;
      index        = NULL;
      flags        = 0;
      Warning      = 0;
      loopNo       = 0;
      tagNo        = 0;
      WrongCat     = NULL;
      WrongTag     = NULL;
      nWrongFields = 0;
    }

    void Data::FreeMemory ( int key )  {
    int i;
      if (name)  delete[] name;
      name = NULL;
      if (Category)  {
        for (i=0;i<nCategories;i++)
          if (Category[i]) delete Category[i];
        delete[] Category;
        Category = NULL;
      }
      nCategories = 0;
      FreeVectorMemory ( index,0 );
      if (key==0)  FreeWrongFields();
    }

    void Data::FreeWrongFields()  {
    int i;
      if (WrongCat)  {
        for (i=0;i<nWrongFields;i++)
          if (WrongCat[i])  delete[] WrongCat[i];
        delete[] WrongCat;
      }
      if (WrongTag)  {
        for (i=0;i<nWrongFields;i++)
          if (WrongTag[i])  delete[] WrongTag[i];
        delete[] WrongTag;
      }
      WrongCat     = NULL;
      WrongTag     = NULL;
      nWrongFields = 0;
    }


    void  Data::SetPrintWarnings ( bool SPW )  {
      if (SPW)  SetFlag    ( CIFFL_PrintWarnings );
          else  RemoveFlag ( CIFFL_PrintWarnings );
    }

    void  Data::SetStopOnWarning ( bool SOW )  {
      if (SOW)  SetFlag    ( CIFFL_StopOnWarnings );
          else  RemoveFlag ( CIFFL_StopOnWarnings );
    }

    void  Data::SetFlag ( CIF_FLAG F )  {
      flags |= F;
    }

    void  Data::RemoveFlag ( CIF_FLAG F )  {
      flags &= ~F;
    }

    void  Data::SetWrongFields ( cpstr *cats, cpstr *tags )  {
    int i,lc,lt;
      FreeWrongFields();
      if ((!cats) || (!tags))  return;
      lc = 0;
      while (cats[lc]) lc++;
      lt = 0;
      while (tags[lt]) lt++;
      nWrongFields = IMax(lc,lt);
      if (nWrongFields>0)  {
        WrongCat = new pstr[nWrongFields];
        WrongTag = new pstr[nWrongFields];
        for (i=0;i<nWrongFields;i++)  {
          WrongCat[i] = NULL;
          WrongTag[i] = NULL;
          if (cats[i])  {
            if (cats[i][0])  CreateCopy ( WrongCat[i],cats[i] );
          }
          if (!WrongCat[i])  {
            CreateCopy ( WrongCat[i],pstr(" ") );
            WrongCat[i][0] = char(1);
          }
          if (tags[i])  CreateCopy ( WrongTag[i],tags[i]  );
                  else  CreateCopy ( WrongTag[i],pstr("") );
        }
      }
    }

    bool Data::CheckWrongField ( cpstr C, cpstr T )  {
    int i;
      for (i=0;i<nWrongFields;i++)
        if ((!strcasecmp(C,WrongCat[i])) &&
            (!strcasecmp(T,WrongTag[i])))  return true;
      return false;
    }

    #define _max_buf_len   500

    static char  _err_string[_max_buf_len+1];
    static int   _err_line;


    int  Data::ReadMMCIFData ( cpstr FName, io::GZ_MODE gzipMode )  {
    io::File f;
    char     S[_max_buf_len+1];
    int      RC,lcount;
      f.assign ( FName,true,false,gzipMode );
      if (f.reset(true))  {
        S[0]   = char(0);
        lcount = 0;
        RC     = ReadMMCIFData ( f,S,lcount );
        f.shut();
        return RC;
      } else  {
        _err_string[0] = char(0);
        _err_line      = 0;
        Warning = CIFRC_CantOpenFile;
        return CIFRC_CantOpenFile;
      }
    }


    // ---------------  General I/O functions

    int  Data::ReadMMCIFData ( io::RFile f, pstr S, int & lcount )  {
    pstr p;
    int  i,llen;
    pstr L;

      FreeMemory(1);
      Warning = 0;
      loopNo  = 0;
      tagNo   = 0;

      // 1. Look for 'data_' tag
      do {
        p = &(S[0]);
        while ((*p==' ') || (*p==char(9)))  p++;
        if (strncmp(p,"data_",5))  {
          f.ReadLine ( S,_max_buf_len );
          lcount++;
          p = NULL;
        }
      } while ((!p) && (!f.FileEnd()));

      if (!p)  {
        strcpy ( _err_string,S );
        _err_line = lcount;
        if (flags & CIFFL_PrintWarnings)
          printf ( "\n **** mmCIF READ ERROR "
                   "<<line %i: no 'data_XXXX' tag found>>\n",lcount );
        return CIFRC_NoDataLine;
      }

      llen = _max_buf_len;
      L    = new char[llen];
      i    = 0;
      p   += 5;
      while ((*p) && (*p!=' ') && (*p!=char(9)))  {
        L[i++] = *p;
        p++;
      }
      L[i] = char(0);
      CreateCopy ( name,L );


      // 2. Loop over tags until next 'data_' or end of file

      while (p)  {

        // skip spaces
        while ((*p==' ') || (*p==char(9)))  p++;

        if ((*p) && (*p!='#'))  {  // this is not a comment, continue
          if (*p=='_')
            GetDataItem ( f,S,L,p,lcount,llen );
          else if (!strncmp(p,"loop_",5))
            GetLoop ( f,S,L,p,lcount,llen );
          else if (!strncmp(p,"data_",5))  {
            p = NULL;
            break;
          } else  {
            // if got to here, the file is corrupted
            strcpy ( _err_string,S );
            _err_line = lcount;
            Warning |= CIFW_UnrecognizedItems;
            if (flags & CIFFL_PrintWarnings)
              printf ( "\n **** mmCIF READ WARNING "
                       "<<line %i: unrecognized string>>\n%s\n",
                       lcount,S );
            while ((*p) && (*p!=' ') && (*p!=char(9)))
              if (*p=='#')  *p = char(0);
                      else  p++;
          }
        } else
          *p = char(0);

        if (Warning && (flags & CIFFL_StopOnWarnings))  {
          if (L)  delete[] L;
          return Warning;
        }

        if (!(*p))  {
          if (!f.FileEnd())  {
            f.ReadLine ( S,_max_buf_len );
            lcount++;
            p = &(S[0]);
          } else
            p = NULL;
        }

      }

      if (L)  delete[] L;

      Optimize();  // get rid of over-allocated fields.

      return Warning;

    }


    void Data::GetDataItem ( io::RFile f, pstr S, pstr & L,
                             pstr & p, int & lcount, int & llen )  {
    PStruct cifStruct;
    char    T[100];
    int     RC,i;

      i = 0;
      while ((*p) && (*p!=' ') && (*p!=char(9)) && (*p!='.'))  {
        if (i<(int)sizeof(T)-1)  T[i++] = *p;
        p++;
      }
      T[i] = char(0);

      if (*p!='.')  {    // category name missing
        strcpy ( L,T );  // item name
        T[0] = char(1);  // special
        T[1] = char(0);  //   category name
      }

      //  look for category
      i = AddCategory ( T );
      if (i<0)  {
        // negative value means that the category was not in the list,
        // but a place for it has been provided and index updated
        cifStruct = new Struct ( T );
        Category[nCategories-1] = cifStruct;
      } else  {
        cifStruct = PStruct(Category[i]);
        if (cifStruct->GetCategoryID()!=MMCIF_Struct)  {
          strcpy ( _err_string,S );
          _err_line = lcount;
          Warning |= CIFW_NotAStructure;
          if (flags & CIFFL_PrintWarnings)
            printf ( "\n **** mmCIF READ WARNING "
                     "<<line %i: %s was a loop -- replaced>>\n%s\n",
                     lcount,T,S );
          delete Category[i];
          cifStruct = new Struct ( T );
          Category[i] = cifStruct;
        }
      }

      if (*p=='.')  {  // get item name
        i = 0;
        p++;  // skip period
        while ((*p) && (*p!=' ') && (*p!=char(9)))  {
          T[i++] = *p;
          p++;
        }
        T[i] = char(0);
      } else
        strcpy ( T,L );

      if (nWrongFields>0)  {
        if (CheckWrongField(cifStruct->name,T))  {
          GetField ( f,S,L,p,lcount,llen );
          cifStruct->DeleteField ( T );
          return;
        }
      }

      RC = GetField ( f,S,L,p,lcount,llen );

      if (RC)  {
        strcpy ( _err_string,S );
        _err_line = lcount;
        Warning |= CIFW_MissingField;
        if (flags & CIFFL_PrintWarnings)
          printf ( "\n **** mmCIF READ WARNING "
                   "<<line %i: expected data field missing>>\n%s\n",
                   lcount,S );
      }

      while ((*p==' ') || (*p==char(9)))  p++;
      if (*p=='#')  *p = char(0);

      i = cifStruct->GetTagNo ( T );
      if (i>=0)  {
        if (flags & CIFFL_SuggestTags)  {
          tagNo++;
          ParamStr ( T,pstr("\1"),tagNo );
        } else  {
          strcpy ( _err_string,S );
          _err_line = lcount;
          Warning |= CIFW_DuplicateTag;
          if (flags & CIFFL_PrintWarnings)
            printf ( "\n **** mmCIF READ WARNING "
                     "<<line %i: duplicated tag>>\n%s\n",lcount,S );
        }
      }
      cifStruct->AddField ( L,T );

    }

    void Data::GetLoop ( io::RFile f, pstr S, pstr & L,
                         pstr & p, int & lcount, int & llen )  {
    PLoop cifLoop;
    pstr  p1;
    char  T[100];
    bool  Repeat,WrongField;
    int   RC,i,nC;

      RC = 0;

      p += 5;  // skip 'loop_' tag

      loopNo++;
      cifLoop = NULL;
      nC      = -1; // undefined category number
      do {

        while ((*p==' ') || (*p==char(9)))  p++;
        p1 = p;

        if (*p=='_')  {

          // get category name
          i = 0;
          while ((*p) && (*p!=' ') && (*p!=char(9)) && (*p!='.'))  {
            if (i<(int)sizeof(T)-1)  T[i++] = *p;
            p++;
          }
          T[i] = char(0);

          if (*p!='.')  {    // category name missing
            strcpy ( L,T );  // item name
            if (flags & CIFFL_SuggestCategories)
                 sprintf ( T,"X%i",loopNo );
            else strcpy  ( T,"X" );
            T[0] = char(1);  // special category name
          }

          if (cifLoop)  {
            if (strcmp(cifLoop->GetCategoryName(),T))  {
              // loop ended, empty
              p    = p1;   // play back to last category
              cifLoop = NULL;
            }
          } else  {
            //  look for category
            i = AddCategory ( T );
            if ((i!=nC) && (nC>=0))  {  // empty loop; most probably
                                        // a corrupted file
              p = p1;   // play back to last category
              strcpy ( _err_string,S );
              _err_line = lcount;
              Warning |= CIFW_EmptyLoop;
              if (flags & CIFFL_PrintWarnings)
                printf ( "\n **** mmCIF READ WARNING "
                         "<<line %i: empty loop>>\n%s\n",lcount,S );
              // AddCategory(..) has added a NULL-Category on the top of
              // category list; remove it now
              DeleteCategory ( nCategories-1 );
              cifLoop = NULL;
    //          return;
            }
            if (i<0)  {
              // negative value means that the category was not in the list,
              // but a place for it has been provided and index updated
              cifLoop = new Loop ( T );
              Category[nCategories-1] = cifLoop;
              nC = nCategories-1;
            }
          }
    /*
     else if (Loop)  {
            if (!strcmp(Loop->GetCategoryName(),
                        Category[i]->GetCategoryName()))  {
              if (Loop->GetCategoryID()!=MMCIF_Loop)  {
                Warning |= CIFW_NotALoop;
                if (flags & CIFFL_PrintWarnings)
                  printf ( "\n **** mmCIF READ WARNING "
                           "<<line %i: %s was a structure --"
                           " replaced>>\n%s\n",lcount,T,S );
                delete Category[i];
                Loop = new Loop ( T );
                Category[i] = Loop;
              }
            } else
              Loop = NULL;
          }
    */
          if (cifLoop)  {

            if (*p=='.')  {  // get item name
              i = 0;
              p++;  // skip period
              while ((*p) && (*p!=' ') && (*p!=char(9)))  {
                T[i++] = *p;
                p++;
              }
              T[i] = char(0);
            } else
              strcpy ( T,L );

            if (nWrongFields>0)
                  WrongField = CheckWrongField ( cifLoop->name,T );
            else  WrongField = false;

            if (!WrongField)  {
              if (cifLoop->AddTag(T)>=0)  {
                if (flags & CIFFL_SuggestTags)  {
                  tagNo++;
                  ParamStr ( T,pstr("\1"),tagNo );
                  cifLoop->AddTag(T);
                } else  {
                  strcpy ( _err_string,S );
                  _err_line = lcount;
                  Warning |= CIFW_DuplicateTag;
                  if (flags & CIFFL_PrintWarnings)
                    printf ( "\n **** mmCIF READ WARNING "
                             "<<line %i: duplicate tag>>\n%s\n",lcount,S );
                }
              }
            }
            Repeat = true;

          } else  {

            p = p1;
            Repeat = false;

          }

        } else if (!(*p) || (*p=='#'))  {
          Repeat = !f.FileEnd();
          if (Repeat)  {
            f.ReadLine ( S,_max_buf_len );
            lcount++;
            p = &(S[0]);
          } else  {
            strcpy ( _err_string,S );
            _err_line = lcount;
            Warning |= CIFW_UnexpectedEOF;
            if (flags & CIFFL_PrintWarnings)
              printf ( "\n **** mmCIF READ WARNING "
                       "<<line %i: unexpected end of file>>\n%s\n",
                       lcount,S );
          }
        } else
          Repeat = false;

      } while (Repeat);

      if (cifLoop)  {
        do  {
          while ((*p==' ') || (*p==char(9)))  p++;
          if (!(*p) || (*p=='#'))  {
            Repeat = !f.FileEnd();
            if (Repeat)  {
              f.ReadLine ( S,_max_buf_len );
              lcount++;
              p = &(S[0]);
            }
          } else if (*p=='_')  Repeat = false;
          else if (!strncmp(p,"loop_",5))  Repeat = false;
          else if (!strncmp(p,"data_",5))  Repeat = false;
          else if (!strncmp(p,"stop_",5))  {
            p += 5;
            Repeat = false;
          } else  {
            RC = GetField ( f,S,L,p,lcount,llen );
            if (!RC)  {
              cifLoop->AddString ( L );
              Repeat = true;
            } else
              Repeat = false;
          }
        } while (Repeat);
        if ((cifLoop->iColumn!=0) || (RC))  {
          strcpy ( _err_string,S );
          _err_line = lcount;
          Warning |= CIFW_LoopFieldMissing;
          if (flags & CIFFL_PrintWarnings)
            printf ( "\n **** mmCIF READ WARNING "
                     "<<line %i: expected loop field missing>>\n%s\n",
                     lcount,S );
        }
      }

    }

    int  Data::GetField ( io::RFile f, pstr S, pstr & L,
                          pstr & p, int & lcount, int & llen )  {
    bool Repeat;
    pstr L1;
    int  i,flen;
    char c;

      flen    = 0;
      L[flen] = char(0);

      do {

        // skip all spaces before the field
        while ((*p==' ') || (*p==char(9)))  p++;

        if ((*p=='#') || (!(*p)))  {
          // comment or end of line met;  the field should be
          // found on the next line
          Repeat = !f.FileEnd();
          if (Repeat)  {
            // take the next line
            f.ReadLine ( S,_max_buf_len );
            lcount++;
            p = &(S[0]);
            Repeat = (*p!=';');
          } else  {
            // end of file and the field is not taken
            L[0] = char(0);
            return 1;
          }
        } else
          // first symbol of a field is found
          Repeat = false;

      } while (Repeat);


      if (*p==';')  {      // this is a multiline field
        p++;
        strcpy ( L,p );    // take first line of the field
        flen = strlen(L);
        while (!f.FileEnd())  {
          f.ReadLine ( S,_max_buf_len );
          lcount++;
          p = &(S[0]);
          if (*p==';')  {  // multiline field terminated
            p++;
            while ((*p==' ') || (*p==char(9)))  p++;
            return 0;
          } else  {        // multiline field continues -- take next line
            flen += strlen(S)+2;
            if (flen>=llen)  {
              llen = flen + IMin(2000,llen);
              L1   = new char[llen];
              strcpy ( L1,L );
              delete[] L;
              L = L1;
            }
            strcat ( L,"\n" );
            strcat ( L,S );
          }
        }

        // end of file -- take the last line of the multiline field
        p = &(S[strlen(S)]);

      } else  {

        i = 0;
        if (*p!='_')  {
          if ((*p=='\'') || (*p=='"'))  {
            c = *p;
            // field in quotation marks
            do  {
              p++;
              // stop taking characters either on end of line
              // or the quotation mark
              while ((*p) && (*p!=c))  {
                L[i++] = *p;
                p++;
              }
              while (*p==c)  {
                // it was a quotation mark -- check that it is followed
                // by end of line or space
                p++;
                if ((*p) && (*p!=' ') && (*p!=char(9)))  {
                  // the quotation mark is not a field terminator and
                  // belongs to the field.
                  L[i++] = c;    // take the quotation mark
                  if (*p!=c)     // take the non-space character
                    L[i++] = *p; // but check for field terminator
                }
              }
              // terminate the loop only on space or end of line
            } while ((*p) && (*p!=' ') && (*p!=char(9)));
            if (*p)  p++;
            L[i] = char(0);
          } else  {
            // a simplest field without spaces
            while ((*p) && (*p!=' ') && (*p!=char(9)))  {
              L[i++] = *p;
              p++;
            }
            L[i] = char(0);
            if (((L[0]=='.') || (L[0]=='?')) && (!L[1]))  {
              //  "no data" tokens
              L[1] = L[0];
              L[0] = char(2);
              L[2] = char(0);
            }
          }
        }

      }

      return 0;

    }

    void  Data::Sort()  {
    int      i,k;
    psvector cnames;

      k = 0;
      for (i=0;i<nCategories;i++)
        if (Category[i])  {
          if (k<i)  Category[k] = Category[i];
          k++;
        }
      for (i=k;i<nCategories;i++)
        Category[i] = NULL;
      nCategories = k;

      FreeVectorMemory ( index ,0 );
      GetVectorMemory  ( cnames,nCategories,0 );
      GetVectorMemory  ( index ,nCategories,0 );

      for (i=0;i<nCategories;i++)  {
        Category[i]->Sort();
        cnames[i] = NULL;
        CreateCopy ( cnames[i],Category[i]->name );
      }

      SortTags ( cnames,nCategories,index );

      for (i=0;i<nCategories;i++)
        if (cnames[i])  delete[] cnames[i];

      if (cnames) delete[] cnames;

    }

    int  Data::GetCategoryNo ( cpstr cname )  {
    //   Binary search for index of category cname in Category[].
    // Return:
    //    >=0 : position of the category found
    //     <0 : the category was not found, it could be inserted before
    //          (-RC-1)th element, where RC is the return value
    int l1,l2,l,k;

      if ((!Category) || (nCategories<1)) return -1;

      if (!index)    Sort();

      if (cname[0])  {
        l  = 0;
        l1 = 0;
        l2 = nCategories-1;
        k  = 1;
        while (l1<l2-1)  {
          l = (l1+l2)/2;
          k = strcasecmp ( cname,Category[index[l]]->name );
          if (k<0)      l2 = l;
          else if (k>0) l1 = l;
          else  {
            l1 = l;
            break;
          }
        }
        if (k==0)  return index[l];    // is at RCth position
        k = strcasecmp(cname,Category[index[l1]]->name);
        if (k==0)  return index[l1];   // is at RCth position
        if (k<0)   return -1;          // would be at (-RC-1)th position
        if (l2!=l1)  {
          k = strcasecmp(cname,Category[index[l2]]->name);
          if (k==0)  return index[l2]; // is at RCth position
          if (k>0)   return -2-l2;   // would be at l2+1=(-RC-1)th position
        }
        return -2-l1;                // would be at l1+1=(-RC-1)th position
      } else
        // 'root' category should be always on top
        if (Category[index[0]]->name[0]==char(1))  return index[0];

      return -1;

    }

    int  Data::AddCategory ( cpstr cname )  {
    //  return -1: a place for category has been added on the top of array;
    //             index is added and sorted automatically
    //        >=0: the category is already in the array -- its position
    //             is returned
    int              l1,l;
    PPCategory Category1;
    ivector          index1;


      if (!Category)  {
        Category     = new PCategory[1];
        Category[0]  = NULL;
        GetVectorMemory ( index,1,0 );
        index[0]     = 0;
        nCategories  = 1;
        return -nCategories;  // the category has been added on the top of array
      }
      l1 = GetCategoryNo ( cname );
      if (l1>=0)  return l1;  // non-negative returns mean that
                              // the category is already in the array
      l1 = -l1-1;  // otherwise the category has to be added and indexed at here
      // put new NULL-category on the top of array and update the index
      Category1 = new PCategory[nCategories+1];
      GetVectorMemory ( index1,nCategories+1,0 );
      for (l=0;l<nCategories;l++)
        Category1[l] = Category[l];
      Category1[nCategories] = NULL;
      for (l=0;l<l1;l++)
        index1[l] = index[l];
      index1[l1] = nCategories;
      for (l=l1+1;l<=nCategories;l++)
        index1[l] = index[l-1];
      delete[] Category;
      FreeVectorMemory ( index,0 );
      Category = Category1;
      index    = index1;
      nCategories++;
      return -nCategories; // the category has been added on
                           // the top of array
    }

    bool Data::WriteMMCIFData ( cpstr FName, io::GZ_MODE gzipMode )  {
    io::File f;
      f.assign ( FName,true,false,gzipMode );
      if (f.rewrite())  {
        WriteMMCIF ( f );
        f.shut();
        return true;
      } else
        return false;
    }

    void Data::WriteMMCIF ( io::RFile f )  {
    int  i;

      if (name)  {
        f.Write     ( pstr("\ndata_") );
        f.WriteLine ( name            );
      } else
        f.WriteLine ( pstr("\ndata_") );

      for (i=0;i<nCategories;i++)
        if (Category[i])
          Category[i]->WriteMMCIF ( f );

    }


    // ---------------  Retrieving data

    int  Data::DeleteCategory ( cpstr CName )  {
    int k;
      k = GetCategoryNo ( CName );
      if (k<0)  return CIFRC_NoCategory;
      return DeleteCategory ( k );
    }

    int  Data::DeleteCategory ( int CatNo ) {
    int i;
      if (Category[CatNo])  delete Category[CatNo];
      for (i=CatNo+1;i<nCategories;i++)
        Category[i-1] = Category[i];
      i = 0;
      while ((i<nCategories) && (index[i]!=CatNo))  {
        if (index[i]>CatNo) index[i]--;
        i++;
      }
      i++;
      while (i<nCategories)  {
        if (index[i]>CatNo) index[i]--;
        index[i-1] = index[i];
        i++;
      }
      nCategories--;
      index   [nCategories] = 0;
      Category[nCategories] = NULL;
      return 0;
    }

    int  Data::DeleteStructure ( cpstr CName )  {
    int k;
      k = GetCategoryNo ( CName );
      if (k<0)  return CIFRC_NoCategory;
      if (Category[k]->GetCategoryID()==MMCIF_Struct)
            return DeleteCategory ( k );
      else  return CIFRC_NotAStructure;
    }

    int  Data::DeleteLoop ( cpstr CName )  {
    int k;
      k = GetCategoryNo ( CName );
      if (k<0)  return CIFRC_NoCategory;
      if (Category[k]->GetCategoryID()==MMCIF_Loop)
            return DeleteCategory ( k );
      else  return CIFRC_NotALoop;
    }

    PCategory Data::GetCategory ( int categoryNo )  {
      if ((categoryNo>=0) && (categoryNo<nCategories))
        return Category[categoryNo];
      return NULL;
    }

    PStruct Data::GetStructure ( cpstr CName )  {
    int i;
      i = GetCategoryNo ( CName );
      if (i<0)  return NULL;
      if (Category[i]->GetCategoryID()==MMCIF_Struct)
            return PStruct(Category[i]);
      else  return NULL;
    }

    PLoop Data::GetLoop ( cpstr CName )  {
    int i;
      i = GetCategoryNo ( CName );
      if (i<0)  return NULL;
      if (Category[i]->GetCategoryID()==MMCIF_Loop)
            return PLoop(Category[i]);
      else  return NULL;
    }

    PLoop Data::FindLoop ( cpstr * tagList )  {
    int i;
      for (i=0;i<nCategories;i++)
        if (Category[i])  {
          if (Category[i]->GetCategoryID()==MMCIF_Loop)  {
            if (Category[i]->CheckTags(tagList))
              return PLoop(Category[i]);
          }
        }
      return NULL;
    }

    void Data::GetDataName ( pstr & dname, bool Remove )  {
      if (Remove)  {
        if (dname)  delete[] dname;
        dname = name;
        name  = NULL;
      } else
        CreateCopy ( dname,name );
    }

    int  Data::CheckData ( cpstr CName, cpstr TName )  {
    //   CheckData(..) returns positive value if the field is in the
    // file:
    //   CIFRC_Structure  category CName is a structure
    //   CIFRC_Loop       category CName is a loop
    // Negative returns mean:
    //   CIFRC_StructureNoTag  category CName is present,
    //                        it is a structure, but it does not
    //                        have tag TName
    //   CIFRC_LoopNoTag       category CName is present,
    //                        it is a loop, but it does not have
    //                        tag TName
    //   CIFRC_NoCategory      category CName is not present.
    // If TName is set to NULL then only the CName is checked and
    // possible returns are CIFRC_Structure, CIFRC_Loop and
    // CIFRC_NoCategory.
    int i,k;
      i = GetCategoryNo ( CName );
      if (i<0)  return CIFRC_NoCategory;
      if (Category[i]->GetCategoryID()==MMCIF_Struct)
            k = CIFRC_Structure;
      else  k = CIFRC_Loop;
      if (TName)  {
        if (Category[i]->GetTagNo(TName)<0)  {
          if (k==CIFRC_Structure)
                k = CIFRC_StructureNoTag;
          else  k = CIFRC_LoopNoTag;
        }
      }
      return k;
    }

    void Data::Optimize()  {
    int              i,k;
    PPCategory C1;
      k = 0;
      for (i=0;i<nCategories;i++)
        if (Category[i])  {
          Category[i]->Optimize();
          if (Category[i]->nTags<=0)  {
            delete Category[i];
            Category[i] = NULL;
          } else
            k++;
        }
      if (k>0)  {
        if (k!=nCategories)  {
          C1 = new PCategory[k];
          k  = 0;
          for (i=0;i<nCategories;i++)
            if (Category[i])
              C1[k++] = Category[i];
          if (Category) delete[] Category;
          Category    = C1;
          nCategories = k;
          FreeVectorMemory ( index,0 );
          Sort();
        }
      } else  {
        if (Category)  delete[] Category;
        Category    = NULL;
        nCategories = 0;
      }
    }

    int  Data::GetString  ( pstr & Dest, cpstr CName,
                                  cpstr TName, bool Remove )  {
    //   GetString(..), GetReal(..) and GetInteger(..) return 0 if the
    // requested field was found and successfully converted. Negative
    // returns mean:
    //    CIFRC_WrongFormat    the field was found but failed to convert
    //                        due to improper numeric format
    //    CIFRC_NoTag          category CName was found, but it does not
    //                        have tag TName
    //    CIFRC_NoCategory     category CName was not found
    //    CIFRC_NotAStructure  category CName was found, but it is a loop
    //                        rather than a structure.
    //   GetString(..) will try to dispose Dest unless it is assugned
    // NULL value before the call. The string will be then dynamically
    // allocated and copied.
    int i = GetCategoryNo ( CName );
      if (i<0)  return CIFRC_NoCategory;
      if (Category[i]->GetCategoryID()!=MMCIF_Struct)
                return CIFRC_NotAStructure;
      return PStruct(Category[i])->GetString ( Dest,TName,Remove );
    }

    pstr Data::GetString ( cpstr CName, cpstr TName, int & RC )  {
    int i = GetCategoryNo ( CName );
      if (i<0)  {
        RC = CIFRC_NoCategory;
        return NULL;
      }
      if (Category[i]->GetCategoryID()!=MMCIF_Struct)  {
        RC = CIFRC_NotAStructure;
        return NULL;
      }
      return PStruct(Category[i])->GetString ( TName,RC );
    }

    int  Data::DeleteField ( cpstr CName, cpstr TName )  {
    int i = GetCategoryNo ( CName );
      if (i<0)  return CIFRC_NoCategory;
      if (Category[i]->GetCategoryID()!=MMCIF_Struct)
                return CIFRC_NotAStructure;
      return PStruct(Category[i])->DeleteField ( TName );
    }

    int  Data::GetReal ( realtype & R, cpstr CName,
                         cpstr TName, bool Remove )  {
    int i = GetCategoryNo ( CName );
      if (i<0)  return CIFRC_NoCategory;
      if (Category[i]->GetCategoryID()!=MMCIF_Struct)
                return CIFRC_NotAStructure;
      return PStruct(Category[i])->GetReal ( R,TName,Remove );
    }

    int  Data::GetInteger ( int & I, cpstr CName,
                                  cpstr TName, bool Remove )  {
    int j = GetCategoryNo ( CName );
      if (j<0)  return CIFRC_NoCategory;
      if (Category[j]->GetCategoryID()!=MMCIF_Struct)
                return CIFRC_NotAStructure;
      return PStruct(Category[j])->GetInteger ( I,TName,Remove );
    }

    int  Data::GetLoopLength ( cpstr CName )  {
    //   GetLoopLength(..) returns CIFRC_NotALoop if the category CName
    // is not a loop, CIFRC_NoCategory if the category CName is not
    // found. Non-negative returns give the length of the loop (may be
    // 0 if the loop is empty).
    int i;
      i = GetCategoryNo ( CName );
      if (i<0)  return CIFRC_NoCategory;
      if (Category[i]->GetCategoryID()!=MMCIF_Loop)
                return CIFRC_NotALoop;
      return PLoop(Category[i])->nRows;
    }

    int  Data::GetLoopString  ( pstr & Dest, cpstr CName,
                                      cpstr TName, int nrow,
                                      bool Remove )  {
    //   GetLoopString(..), GetLoopReal(..) and GetLoopInteger(..) act
    // like GetString(..), GetReal(..) and GetInteger(..) above for
    // nrow-th element of the 'loop_' (indexed like 0..N-1 where N
    // is obtained through GetLoopLength(..)). They will return
    // CIFRC_WrongIndex if nrow is out of range.
    int i = GetCategoryNo ( CName );
      if (i<0)  return CIFRC_NoCategory;
      if (Category[i]->GetCategoryID()!=MMCIF_Loop)
                return CIFRC_NotALoop;
      return PLoop(Category[i])->GetString ( Dest,TName,nrow,Remove );
    }

    pstr Data::GetLoopString  ( cpstr CName, cpstr TName,
                                      int nrow, int & RC )  {
    int i = GetCategoryNo ( CName );
      if (i<0)  {
        RC = CIFRC_NoCategory;
        return NULL;
      }
      if (Category[i]->GetCategoryID()!=MMCIF_Loop)  {
        RC = CIFRC_NotALoop;
        return NULL;
      }
      return  PLoop(Category[i])->GetString ( TName,nrow,RC );
    }

    int Data::DeleteLoopField ( cpstr CName, cpstr TName,
                                      int nrow )  {
    int i = GetCategoryNo ( CName );
      if (i<0)  return CIFRC_NoCategory;
      if (Category[i]->GetCategoryID()!=MMCIF_Loop)
                return CIFRC_NotALoop;
      return PLoop(Category[i])->DeleteField ( TName,nrow );
    }

    int  Data::GetLoopReal ( realtype & R, cpstr CName,
                                   cpstr TName, int nrow,
                                   bool Remove )  {
    int i = GetCategoryNo ( CName );
      if (i<0)  return CIFRC_NoCategory;
      if (Category[i]->GetCategoryID()!=MMCIF_Loop)
                return CIFRC_NotALoop;
      return PLoop(Category[i])->GetReal ( R,TName,nrow,Remove );
    }

    int  Data::GetLoopInteger ( int & I, cpstr CName,
                                      cpstr TName, int nrow,
                                      bool Remove )  {
    int j = GetCategoryNo ( CName );
      if (j<0)  return CIFRC_NoCategory;
      if (Category[j]->GetCategoryID()!=MMCIF_Loop)
                return CIFRC_NotALoop;
      return PLoop(Category[j])->GetInteger ( I,TName,nrow,Remove );
    }


    int  Data::GetLoopSVector ( psvector & S, cpstr CName,
                                      cpstr TName, int i1, int i2,
                                      bool Remove )  {
    //   GetLoopSVector(..), GetLoopRVector(..) and GetLoopIVector(..)
    // read CIF 'loop_' data into allocated vectors of strings, reals and
    // integers, correspondingly. The vectors may be deallocated prior
    // to call and assigned NULL, in which case they will be allocated
    // with offsets of i1, which is also the lower index of the 'loop_'
    // data transferred into it. The upper vector index is given by i2 or
    // by the loop's length whichever is less. If vectors are not assigned
    // NULL prior the call, it is assumed that they are properly (i1-offset,
    // i2-i1+1 length) allocated.
    //   The return codes are same as those of GetLoopString(..),
    // GetLoopReal(..) and GetLoopInteger(..).
    int i;
      i = GetCategoryNo ( CName );
      if (i<0)    return CIFRC_NoCategory;
      if (Category[i]->GetCategoryID()!=MMCIF_Loop)
                  return CIFRC_NotALoop;
      return PLoop(Category[i])->GetSVector ( S,TName,i1,i2,Remove );
    }

    int  Data::GetLoopRVector ( rvector & R, cpstr CName,
                                      cpstr TName, int i1, int i2,
                                      bool Remove )  {
    int i;
      i = GetCategoryNo ( CName );
      if (i<0)    return CIFRC_NoCategory;
      if (Category[i]->GetCategoryID()!=MMCIF_Loop)
                  return CIFRC_NotALoop;
      return PLoop(Category[i])->GetRVector ( R,TName,i1,i2,Remove );
    }

    int  Data::GetLoopIVector ( ivector & I, cpstr CName,
                                      cpstr TName, int i1, int i2,
                                      bool Remove )  {
    int j;
      j = GetCategoryNo ( CName );
      if (j<0)    return CIFRC_NoCategory;
      if (Category[j]->GetCategoryID()!=MMCIF_Loop)
                  return CIFRC_NotALoop;
      return PLoop(Category[j])->GetIVector ( I,TName,i1,i2,Remove );
    }


    // ---------------  Storing data

    void  Data::PutDataName ( cpstr dname )  {
    // stores name for 'data_' record
      CreateCopy ( name,dname );
    }

    int  Data::PutNoData ( int NoDataType, cpstr CName, cpstr TName )  {
    PStruct cifStruct;
    int     i,RC;

      RC = CIFRC_Ok;

      //  look for category
      i = AddCategory ( CName );
      if (i<0)  {
        // negative value means that the category was not in the list,
        // but a place for it has been provided and index updated
        cifStruct = new Struct ( CName );
        Category[nCategories-1] = cifStruct;
      } else  {
        cifStruct = PStruct(Category[i]);
        if (cifStruct->GetCategoryID()!=MMCIF_Struct)  {
          RC = CIFRC_NotAStructure;
          delete Category[i];
          cifStruct = new Struct ( CName );
          Category[i] = cifStruct;
        }
      }

      cifStruct->PutNoData ( NoDataType,TName );

      return RC;

    }

    int  Data::PutString ( cpstr S, cpstr CName,
                           cpstr TName, bool Concatenate )  {
    //   PutString(..), PutReal(..) and PutInteger(..) will put the
    // values given into the specified category (CName) under the
    // specified tag (TName). The category, tag and field are created
    // automatically; the field will be replaced silently if identical
    // CName.TName is specified in two calls. Calls of these functions
    // may follow in random order; however CIF file will have all tags
    // grouped by categories and catgories will follow in the order
    // of first appearance in PutString(..), PutReal(..) or
    // PutInteger(..).
    //   Return code - one of CIFRC_Ok or CIFRC_NotAStruct
    PStruct cifStruct;
    int     i,RC;

      RC = CIFRC_Ok;

      //  look for category
      i = AddCategory ( CName );
      if (i<0)  {
        // negative value means that the category was not in the list,
        // but a place for it has been provided and index updated
        cifStruct = new Struct ( CName );
        Category[nCategories-1] = cifStruct;
      } else  {
        cifStruct = PStruct(Category[i]);
        if (cifStruct->GetCategoryID()!=MMCIF_Struct)  {
          RC = CIFRC_NotAStructure;
          delete Category[i];
          cifStruct = new Struct ( CName );
          Category[i] = cifStruct;
        }
      }

      cifStruct->AddField ( S,TName,Concatenate );

      return RC;

    }

    int   Data::PutDate ( cpstr CName, cpstr TName )  {
    PStruct cifStruct;
    int     i,RC;

      RC = CIFRC_Ok;

      //  look for category
      i = AddCategory ( CName );
      if (i<0)  {
        // negative value means that the category was not in the list,
        // but a place for it has been provided and index updated
        cifStruct = new Struct ( CName );
        Category[nCategories-1] = cifStruct;
      } else  {
        cifStruct = PStruct(Category[i]);
        if (cifStruct->GetCategoryID()!=MMCIF_Struct)  {
          RC = CIFRC_NotAStructure;
          delete Category[i];
          cifStruct = new Struct ( CName );
          Category[i] = cifStruct;
        }
      }

      cifStruct->PutDate ( TName );

      return RC;

    }

    int  Data::PutReal ( realtype R, cpstr CName,
                         cpstr TName, int prec )  {
    char rS[100];
      sprintf ( rS,"%.*g",prec,R );
      return PutString ( rS,CName,TName,false );
    }

    int  Data::PutInteger ( int I, cpstr CName,
                                         cpstr TName )  {
    char iS[100];
      if (I>MinInt4)  sprintf ( iS,"%i",I );
      else  {
        iS[0] = char(2);
        iS[1] = '.';
        iS[2] = char(0);
      }
      return PutString ( iS,CName,TName,false );
    }


    int  Data::AddLoop ( cpstr CName, PLoop & cifLoop )  {
    int  i,RC;

      RC = CIFRC_Ok;

      //  look for category
      i = AddCategory ( CName );
      if (i<0)  {
        // negative value means that the category was not in the list,
        // but a place for it has been provided and index updated
        cifLoop = new Loop ( CName );
        Category[nCategories-1] = cifLoop;
        RC = CIFRC_Created;
      } else  {
        cifLoop = PLoop(Category[i]);
        if (cifLoop->GetCategoryID()!=MMCIF_Loop)  {
          RC = CIFRC_NotALoop;
          delete Category[i];
          cifLoop = new Loop ( CName );
          Category[i] = cifLoop;
        }
      }

      return RC;

    }

    int  Data::AddStructure ( cpstr CName, PStruct & cifStruct )  {
    int  i,RC;

      RC = CIFRC_Ok;

      //  look for category
      i = AddCategory ( CName );
      if (i<0)  {
        // negative value means that the category was not in the list,
        // but a place for it has been provided and index updated
        cifStruct = new Struct ( CName );
        Category[nCategories-1] = cifStruct;
        RC = CIFRC_Created;
      } else  {
        cifStruct = PStruct(Category[i]);
        if (cifStruct->GetCategoryID()!=MMCIF_Struct)  {
          RC = CIFRC_NotAStructure;
          delete Category[i];
          cifStruct = new Struct ( CName );
          Category[i] = cifStruct;
        }
      }

      return RC;

    }


    int  Data::PutLoopNoData ( int NoDataType, cpstr CName,
                                     cpstr TName, int nrow )  {
    PLoop cifLoop;
    int   i,RC;

      RC = CIFRC_Ok;

      //  look for category
      i = AddCategory ( CName );
      if (i<0)  {
        // negative value means that the category was not in the list,
        // but a place for it has been provided and index updated
        cifLoop = new Loop ( CName );
        Category[nCategories-1] = cifLoop;
      } else  {
        cifLoop = PLoop(Category[i]);
        if (cifLoop->GetCategoryID()!=MMCIF_Loop)  {
          RC = CIFRC_NotALoop;
          delete Category[i];
          cifLoop = new Loop ( CName );
          Category[i] = cifLoop;
        }
      }

      cifLoop->PutNoData ( NoDataType,TName,nrow );

      return RC;

    }


    int  Data::PutLoopString ( cpstr S, cpstr CName,
                                     cpstr TName, int nrow )  {
    //   PutLoopString(..), PutLoopReal(..) and PutLoopInteger(..) act
    // like PutString(..), PutReal(..) and PutInteger(..) above for
    // nrow-th element of the 'loop_' CName (indexed begining from 0).
    // In consequitive calls, given values of nrow does not have to be
    // ordered; the most efficient way is to start with HIGHEST value
    // for nrow in the loop and move down to 0. The least efficient way
    // is to start with nrow=0 and move up.
    PLoop cifLoop;
    int   i,RC;

      RC = CIFRC_Ok;

      //  look for category
      i = AddCategory ( CName );
      if (i<0)  {
        // negative value means that the category was not in the list,
        // but a place for it has been provided and index updated
        cifLoop = new Loop ( CName );
        Category[nCategories-1] = cifLoop;
      } else  {
        cifLoop = PLoop(Category[i]);
        if (cifLoop->GetCategoryID()!=MMCIF_Loop)  {
          RC = CIFRC_NotALoop;
          delete Category[i];
          cifLoop = new Loop ( CName );
          Category[i] = cifLoop;
        }
      }

      cifLoop->PutString ( S,TName,nrow );

      return RC;

    }

    int  Data::PutLoopReal ( realtype R, cpstr CName,
                                   cpstr TName, int nrow, int prec )  {
    char rS[100];
      sprintf ( rS,"%.*g",prec,R );
      return PutLoopString ( rS,CName,TName,nrow );
    }

    int  Data::PutLoopInteger ( int I, cpstr CName,
                                      cpstr TName, int nrow )  {
    char iS[100];
      sprintf ( iS,"%i",I );
      return PutLoopString ( iS,CName,TName,nrow );
    }


    int  Data::PutLoopSVector ( psvector S, cpstr CName,
                                      cpstr TName, int i1, int i2 )  {
    //   PutLoopSVector(..), PutLoopRVector(..) and PutLoopIVector(..)
    // put vectors of values into specified loop fields. Parameters i1
    // and i2 give the range of indices of values which are to be
    // transfered. To transfer an entire vector allocated as [0..N-1]
    // i1 shoudl be set to 0 and i2 - to N-1. Note that the loop is
    // always indexed as starting form 0 on, therefore negative i1 and
    // i2 are not allowed, and specifying i1>0 will leave first i1
    // elements of the CIF loop for the corresponding tag undefined
    // (will be output like '?').
    PLoop cifLoop;
    int   i,RC;

      RC = CIFRC_Ok;

      //  look for category
      i = AddCategory ( CName );
      if (i<0)  {
        // negative value means that the category was not in the list,
        // but a place for it has been provided and index updated
        cifLoop = new Loop ( CName );
        Category[nCategories-1] = cifLoop;
      } else  {
        cifLoop = PLoop(Category[i]);
        if (cifLoop->GetCategoryID()!=MMCIF_Loop)  {
          RC = CIFRC_NotALoop;
          delete Category[i];
          cifLoop = new Loop ( CName );
          Category[i] = cifLoop;
        }
      }

      cifLoop->PutSVector ( S,TName,i1,i2 );

      return RC;

    }

    int  Data::PutLoopRVector ( rvector R, cpstr CName,
                                      cpstr TName,
                                      int i1, int i2, int prec )  {
    PLoop cifLoop;
    int   i,RC;

      RC = CIFRC_Ok;

      //  look for category
      i = AddCategory ( CName );
      if (i<0)  {
        // negative value means that the category was not in the list,
        // but a place for it has been provided and index updated
        cifLoop = new Loop ( CName );
        Category[nCategories-1] = cifLoop;
      } else  {
        cifLoop = PLoop(Category[i]);
        if (cifLoop->GetCategoryID()!=MMCIF_Loop)  {
          RC = CIFRC_NotALoop;
          delete Category[i];
          cifLoop = new Loop ( CName );
          Category[i] = cifLoop;
        }
      }

      cifLoop->PutRVector ( R,TName,i1,i2,prec );

      return RC;

    }

    int  Data::PutLoopIVector ( ivector I, cpstr CName,
                                      cpstr TName,
                                      int i1, int i2 )  {
    PLoop cifLoop;
    int   j,RC;

      RC = CIFRC_Ok;

      //  look for category
      j = AddCategory ( CName );
      if (j<0)  {
        // negative value means that the category was not in the list,
        // but a place for it has been provided and index updated
        cifLoop = new Loop ( CName );
        Category[nCategories-1] = cifLoop;
      } else  {
        cifLoop = PLoop(Category[j]);
        if (cifLoop->GetCategoryID()!=MMCIF_Loop)  {
          RC = CIFRC_NotALoop;
          delete Category[j];
          cifLoop = new Loop ( CName );
          Category[j] = cifLoop;
        }
      }

      cifLoop->PutIVector ( I,TName,i1,i2 );

      return RC;

    }


    int Data::RenameCategory ( cpstr CName,
                                     cpstr newCName )  {
    int i,RC;
      i = GetCategoryNo ( CName );
      if (i>=0)  {
        Category[i]->PutCategoryName ( newCName );
        Sort();
        RC = CIFRC_Ok;
      } else
        RC = CIFRC_NoCategory;
      return RC;
    }


    void Data::Copy ( PData Data )  {
    int i;
      FreeMemory(0);
      CreateCopy ( name,Data->name );
      nCategories = Data->nCategories;
      if (nCategories>0)  {
        Category = new PCategory[nCategories];
        GetVectorMemory ( index,nCategories,0 );
        for (i=0;i<nCategories;i++)  {
          if (Data->Category[i])  {
            if (Data->Category[i]->GetCategoryID()==MMCIF_Struct)
                  Category[i] = new Struct();
            else  Category[i] = new Loop();
            Category[i]->Copy ( Data->Category[i] );
          } else
            Category[i] = NULL;
          index[i] = Data->index[i];
        }
      }
      flags   = Data->flags;
      Warning = Data->Warning;
    }


    int  Data::CopyCategory ( PData Data, cpstr CName,
                                    cpstr newCName )  {
    PCategory Cat;
    int             i,di,dc,RC;

      di = Data->GetCategoryNo ( CName );

      if (di>=0)  {

        RC = CIFRC_Ok;
        dc = Data->Category[di]->GetCategoryID();

        //  look for category
        i = AddCategory ( CName );
        if (i<0)  {
          // negative value means that the category was not in the list,
          // but a place for it has been provided and index updated
          if (dc==MMCIF_Loop)  Cat = new Loop   ( CName );
                         else  Cat = new Struct ( CName );
          Category[nCategories-1] = Cat;
        } else  {
          Cat = Category[i];
          if (Cat->GetCategoryID()!=dc)  {
            RC = CIFRC_NotAStructure;
            delete Category[i];
            if (dc==MMCIF_Loop)  Cat = new Loop   ( CName );
                           else  Cat = new Struct ( CName );
            Category[i] = Cat;
          }
        }

        Cat->Copy ( Data->Category[di] );
        if (newCName)  {
          Cat->PutCategoryName ( newCName );
          Sort();
        }

      } else
        RC = CIFRC_NoCategory;

      return RC;

    }

    void Data::PrintCategories()  {
    // for debuging only
    int i;
      printf ( " Total %i categories:\n",nCategories );
      for (i=0;i<nCategories;i++)
        if (Category[i])  {
          printf ( " %5i. ",i+1 );
          if (Category[i]->GetCategoryID()==MMCIF_Loop)
                printf ( "Loop      %s\n",Category[i]->name );
          else  printf ( "Structure %s\n",Category[i]->name );
        }
    }


    void Data::write ( io::RFile f )  {
    int i,k;
      if (!index)  Sort();
      f.CreateWrite ( name );
      f.WriteInt    ( &nCategories );
      for (i=0;i<nCategories;i++)  {
        if (Category[i])  {
          k = Category[i]->GetCategoryID();
          f.WriteInt ( &k );
          Category[i]->write ( f );
        } else  {
          k = -1;
          f.WriteInt ( &k );
        }
        f.WriteInt ( &(index[i]) );
      }
      f.WriteInt ( &flags   );
      f.WriteInt ( &Warning );
    }


    void Data::read ( io::RFile f )  {
    int i,k;
      FreeMemory(0);
      f.CreateRead ( name );
      f.ReadInt    ( &nCategories );
      if (nCategories>0)  {
        Category = new PCategory[nCategories];
        GetVectorMemory ( index,nCategories,0 );
        for (i=0;i<nCategories;i++)  {
          f.ReadInt ( &k );
          if (k>=0)  {
            if (k==MMCIF_Struct)  Category[i] = new Struct();
                            else  Category[i] = new Loop();
            Category[i]->read ( f );
          } else
            Category[i] = NULL;
          f.ReadInt ( &(index[i]) );
        }
      }
      f.ReadInt ( &flags   );
      f.ReadInt ( &Warning );
    }


    MakeStreamFunctions(Data)



    //  ======================  File  =============================


    File::File() : io::Stream()  {
      InitFile();
    }

    File::File ( cpstr FName, io::GZ_MODE gzipMode ) : io::Stream()  {
      InitFile ();
      ReadMMCIFFile ( FName,gzipMode );
    }

    File::File ( io::RPStream Object ) : io::Stream(Object)  {
      InitFile();
    }

    File::~File()  {
      FreeMemory();
    }

    void File::InitFile()  {
      nData         = 0;
      nAllocData    = 0;
      data          = NULL;
      index         = NULL;
      PrintWarnings = false;
      StopOnWarning = false;
    }

    void File::FreeMemory()  {
    int i;
      for (i=0;i<nData;i++)
        if (data[i])  delete data[i];
      if (data)  delete[] data;
      data       = NULL;
      FreeVectorMemory ( index,0 );
      nData      = 0;
      nAllocData = 0;
    }


    pstr GetMMCIFInputBuffer ( int & LineNo )  {
      LineNo = _err_line;
      _err_string[sizeof(_err_string)-1] = char(0);
      return _err_string;
    }

    int  File::ReadMMCIFFile ( cpstr FName, io::GZ_MODE gzipMode )  {
    io::File  f;
    PData     CIF;
    char      S[500];
    int       RC,lcount;

      FreeMemory();

      CIF = NULL;
      f.assign ( FName,true,false,gzipMode );
      if (f.reset(true))  {
        S[0]   = char(0);
        lcount = 0;
        RC     = 0;
        while ((!RC) && (!f.FileEnd()))  {
          if (!CIF)  CIF = new Data();
          CIF->SetPrintWarnings ( PrintWarnings );
          CIF->SetStopOnWarning ( StopOnWarning );
          RC = CIF->ReadMMCIFData ( f,S,lcount );
          if (!RC)  {
            ExpandData ( nData+1 );
            data[nData] = CIF;
            nData++;
            CIF = NULL;
          }
        }
        if (CIF)  delete CIF;
        f.shut();
        if (RC==CIFRC_NoDataLine)  {
          if (nData>0)  RC = 0;
        }
        Sort();
        return RC;
      } else
        return CIFRC_CantOpenFile;

    }

    int File::WriteMMCIFFile ( cpstr FName, io::GZ_MODE gzipMode )  {
    io::File f;
      f.assign ( FName,true,false,gzipMode );
      if (f.rewrite())  {
        WriteMMCIF ( f );
        f.shut();
        return 0;
      } else
        return CIFRC_CantOpenFile;
    }

    void File::WriteMMCIF ( io::RFile f )  {
    int i;
      for (i=0;i<nData;i++)
        if (data[i])
          data[i]->WriteMMCIF ( f );
    }


    void File::Sort()  {
    psvector tag;
    int      i;
      if (nData>0)  {
        FreeVectorMemory ( index,0 );
        GetVectorMemory  ( index,nData,0 );
        GetVectorMemory  ( tag  ,nData,0 );
        for (i=0;i<nData;i++)  {
          tag[i] = NULL;
          CreateCopy ( tag[i],data[i]->name );
        }
        SortTags ( tag,nData,index );
        for (i=0;i<nData;i++)
          if (tag[i])  {
            delete[] tag[i];
            tag[i] = NULL;
          }
        FreeVectorMemory ( tag,0 );
      }
    }

    int  File::GetCIFDataNo ( cpstr DName )  {
    //   Binary search for index of DName ttag in data[].
    // Return:
    //    >=0 : position of the DName found
    //     <0 : the DName was not found, it could be inserted before
    //          (-RC-1)th element, where RC is the return value
    int l1,l2,l,k;

      if (!data)   return -1;
      if (!index)  Sort();

      l  = 0;
      l1 = 0;
      l2 = nData-1;
      k  = 1;
      while (l1<l2-1)  {
        l = (l1+l2)/2;
        k = strcasecmp ( DName,data[index[l]]->name );
        if (k<0)      l2 = l;
        else if (k>0) l1 = l;
        else  {
          l1 = l;
          break;
        }
      }

      if (k==0)  return index[l];    // is at RCth position
      k = strcasecmp ( DName,data[index[l1]]->name );
      if (k==0)  return index[l1];   // is at RCth position
      if (k<0)   return -1;          // would be at (-RC-1)th position
      if (l2!=l1)  {
        k = strcasecmp ( DName,data[index[l2]]->name );
        if (k==0)  return index[l2]; // is at RCth position
        if (k>0)   return -2-l2;     // would be at l2+1=(-RC-1)th position
      }

      return -2-l1;                  // would be at l1+1=(-RC-1)th position

    }

    PData  File::GetCIFData ( int dataNo )  {
      if ((dataNo>=0) && (dataNo<nData))  return data[dataNo];
                                    else  return NULL;
    }

    PData  File::GetCIFData ( cpstr DName )  {
    int l;
      l = GetCIFDataNo ( DName );
      if (l>=0)  return data[l];
           else  return NULL;
    }

    void File::ExpandData ( int nDataNew )  {
    int          i,nAD;
    PPData data1;
    ivector      index1;
      if (nDataNew>nAllocData)  {
        nAD   = nDataNew + IMin(nAllocData/2+1,100);
        data1 = new PData[nAD];
        GetVectorMemory ( index1,nAD,0 );
        for (i=0;i<nAllocData;i++)  {
          data1 [i] = data [i];
          index1[i] = index[i];
        }
        for (i=nAllocData;i<nAD;i++)  {
          data1 [i] = NULL;
          index1[i] = i;
        }
        if (data)  delete[] data;
        FreeVectorMemory ( index,0 );
        data       = data1;
        index      = index1;
        nAllocData = nAD;
      }
    }

    int  File::AddCIFData ( cpstr DName )  {
    //  return -1: the CIF data structure has been added on the
    //             top of data array; the index is added and sorted
    //             automatically
    //        >=0: the CIF data structure is already in the array
    //             -- its position is returned
    int  i1,i;
      if (!data)  {
        ExpandData ( 3 );  // get space for first 3 CIF data structures
        data[0] = new Data ( DName );
        nData   = 1;
        return -nData;     // the CIF data structure has been added
                           // "on the top" of array
      }
      i1 = GetCIFDataNo ( DName );
      if (i1>=0)  return i1;  // non-negative returns mean that the CIF
                              // data structure is already in the array
      i1 = -i1-1;  // otherwise the data has to be added and indexed at here
      // put new CIF data structure on the top of array and update index
      ExpandData ( nData+1 );
      data[nData] = new Data ( DName );
      for (i=nData;i>i1;i--)
        index[i] = index[i-1];
      index[i1] = nData;
      nData++;
      return -nData; // the tag has been added on the top of array
    }


    int File::DeleteCIFData ( cpstr DName )  {
    int dataNo = GetCIFDataNo ( DName );

      if (dataNo>=0)  return DeleteCIFData ( dataNo );
      return dataNo;

    }

    int File::DeleteCIFData ( int dataNo )  {
    int i;

      if ((0<=dataNo) && (dataNo<nData))  {

        if (data[dataNo])  delete data[dataNo];
        for (i=dataNo+1;i<nData;i++)
          data[i-1] = data[i];
        nData--;

        Sort();

        return 0;

      }

      return -nData;

    }

    void File::Copy  ( PFile File )  {
    int i;
      FreeMemory();
      nData      = File->nData;
      nAllocData = nData;
      if (nData>0)  {
        data = new PData[nData];
        for (i=0;i<nData;i++)  {
          if (File->data[i])  {
            data[i] = new Data();
            data[i]->Copy ( File->data[i] );
          } else
            data[i] = NULL;
        }
      }
    }


    void File::write ( io::RFile f )  {
    int i,k;
      f.WriteInt ( &nData );
      for (i=0;i<nData;i++)
        if (data[i])  {
          k = 1;
          f.WriteInt ( &k );
          data[i]->write ( f );
        } else  {
          k = 0;
          f.WriteInt ( &k );
        }
    }

    void File::read  ( io::RFile f )  {
    int i,k;
      FreeMemory();
      f.ReadInt ( &nData );
      nAllocData = nData;
      if (nData>0)  {
        data = new PData[nData];
        for (i=0;i<nData;i++)  {
          f.ReadInt ( &k );
          if (k)  {
            data[i] = new Data();
            data[i]->read ( f );
          } else
            data[i] = NULL;
        }
      }
    }


    MakeStreamFunctions(File)


    int  isCIF ( cpstr FName, io::GZ_MODE gzipMode )  {
    io::File f;
    int      rc;

      f.assign ( FName,true,false,gzipMode );
      if (f.reset(true))  {
        rc = isCIF ( f );
        f.shut();
      } else
        rc = -1;

      return rc;

    }

    int  isCIF ( io::RFile f )  {
    char    S[_max_buf_len+1];
    bool Done;
    pstr    p;

      f.ReadLine ( S,_max_buf_len );
      S[_max_buf_len] = char(0);
      Done = false;
      while (!Done)  {
        p = &(S[0]);
        while ((*p==' ') || (*p==char(9)))  p++;
        Done = !strncmp(p,"data_",5);
        if (!Done)  {
          if (f.FileEnd())  {
            Done = true;
            p    = NULL;
          } else  {
            f.ReadLine ( S,_max_buf_len );
            S[_max_buf_len] = char(0);
          }
        }
      }

      if (!p)  return 1;
      if (!strncmp(p,"data_",5))  return 0;
                            else  return 1;

    }


    pstr GetCIFMessage ( pstr M, int RC )  {
    int  LineNo;
    pstr InputLine;

      InputLine = GetMMCIFInputBuffer ( LineNo );

      if (RC>10)  {
        if (RC & CIFW_UnrecognizedItems)
          sprintf ( M,"unrecognized items found on %ith line\n%s",
                    LineNo,InputLine );
        else if (RC & CIFW_MissingField)
          sprintf ( M,"expected data field not found; line %i reads\n%s",
                    LineNo,InputLine );
        else if (RC & CIFW_EmptyLoop)
          sprintf ( M,"empty loop ('loop_') on %ith line\n%s",
                    LineNo,InputLine );
        else if (RC & CIFW_UnexpectedEOF)
          sprintf ( M,"unexpected end of file; line %i reads\n%s",
                    LineNo,InputLine );
        else if (RC & CIFW_LoopFieldMissing)
          sprintf ( M,"expected data field in a loop not found; "
                      "line %i reads\n%s", LineNo,InputLine );
        else if (RC & CIFW_LoopFieldMissing)
          sprintf ( M,"expected data field in a loop not found; "
                      "line %i reads\n%s", LineNo,InputLine );
        else if (RC & CIFW_NotAStructure)
          sprintf ( M,"a loop is used as a structure on line %i\n%s",
                    LineNo,InputLine );
        else if (RC & CIFW_NotALoop)
          sprintf ( M,"a structure is used as a loop on line %i\n%s",
                    LineNo,InputLine );
        else if (RC & CIFW_DuplicateTag)
          sprintf ( M,"duplicate tag was found on line %i\n%s",
                    LineNo,InputLine );
        else
          sprintf ( M,"undocumented warning issued for line %i\n%s",
                    LineNo,InputLine );
      } else if (RC<0)
        switch (RC)  {
          case CIFRC_StructureNoTag : strcpy(M,"tag of a structure not "
                                               "found");
                                 break;
          case CIFRC_LoopNoTag      : strcpy(M,"tag of a loop not found");
                                 break;
          case CIFRC_NoCategory     : strcpy(M,"category not found");
                                 break;
          case CIFRC_WrongFormat    : strcpy(M,"wrong format of a number");
                                 break;
          case CIFRC_NoTag          : strcpy(M,"tag not found");
                                 break;
          case CIFRC_NotAStructure  : strcpy(M,"category is not a "
                                               "structure");
                                 break;
          case CIFRC_NotALoop       : strcpy(M,"category is not a loop");
                                 break;
          case CIFRC_WrongIndex     : strcpy(M,"index outside the loop's "
                                               "limits");
                                 break;
          case CIFRC_NoField        : strcpy(M,"data is absent");
                                 break;
          case CIFRC_Created        : strcpy(M,"category created");
                                 break;
          case CIFRC_CantOpenFile   : strcpy(M,"can't open CIF file");
                                 break;
          case CIFRC_NoDataLine     : strcpy(M,"'data_' tag not found." );
                                 break;
          default                   : strcpy(M,"undocumented return code");
        }

      return M;

    }

  }  // namespace mmcif

}  // namespace mmdb
