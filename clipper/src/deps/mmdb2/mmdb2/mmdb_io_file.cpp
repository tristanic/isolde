//  $Id: file_.cpp,v 1.29 2012/01/26 17:52:19 ekr Exp $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2008.
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
//    03.12.15   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  mmdb_io_file  <implementation>
//       ~~~~~~~~~
//  **** Classes :  mmdb::io::File  - file I/O Support.
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2000-2015
//
//  =================================================================
//


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef  _WIN32
# include <windows.h>
# define sleep Sleep
#endif

#if !defined _WIN32 || defined __MINGW32__
# ifndef  __UNISTD_H
#  include <unistd.h>
# endif
#endif


#include "mmdb_io_file.h"


// _WIN32_NEWLINE should be raised when compilinig on Windows in
// order to enforce Windows' line endings when writing text in
// files opened for *binary* output. Otherwise, writing text lines
// in binary files will results in UNIX line endings. Line endings
// in files, opened for output in text mode, will be always
// platform-specific.

#ifdef  _WIN32_NEWLINE
//  for DOS/WINDOWS machines:
#define  NEWLINE  "\r\n"
#else
//  for UNIX machines:
#define  NEWLINE  "\n"
#endif

namespace mmdb  {

  namespace io  {

#ifdef _WIN32
    const char _dir_sep_c = '\\';
    cpstr      _dir_sep   = "\\";
#else
    const char _dir_sep_c = '/';
    cpstr      _dir_sep   = "/";
#endif

    // ===================  Auxilary Functions  =========================

    cpstr GetFPath ( pstr FilePath, SYSKEY syskey )  {
    pstr P;

      if (syskey==syskey_unix)
        P = LastOccurence(FilePath,'/');
      else if (syskey==syskey_win)
        P = LastOccurence(FilePath,'\\');
      else if (syskey==syskey_all)  {
        P = LastOccurence(FilePath,'/');
        if (!P)  P = LastOccurence(FilePath,'\\');
      } else
        P = NULL;

      if (P)  {
        P = P + 1;
        *P = char(0);
      } else
        FilePath[0] = char(0);

      return FilePath;

    }

    cpstr GetFName ( cpstr FilePath, SYSKEY syskey )  {
    pstr P;
      if (syskey==syskey_unix)
        P = LastOccurence(FilePath,'/');
      else if (syskey==syskey_win)
        P = LastOccurence(FilePath,'\\');
      else if (syskey==syskey_all)  {
        P = LastOccurence(FilePath,'/');
        if (!P)  P = LastOccurence(FilePath,'\\');
      } else
        P = NULL;
      if (!P)  return FilePath;
         else  return P + 1;
    }

    cpstr GetFExt ( cpstr FilePath )  {
    pstr P;
      P = FirstOccurence ( GetFName(FilePath),'.');
      if (!P) return &(FilePath[strlen(FilePath)]);
         else return P;
    }

    cpstr ChangeExt ( pstr FilePath, cpstr newExt, SYSKEY syskey )  {
    int i;
      i = strlen(FilePath)-1;
      if (syskey==syskey_unix)
        while ((i>0) && (FilePath[i]!='.') && (FilePath[i]!='/'))  i--;
      else if (syskey==syskey_win)
        while ((i>0) && (FilePath[i]!='.') && (FilePath[i]!='\\'))  i--;
      else if (syskey==syskey_all)
        while ((i>0) && (FilePath[i]!='.') &&
               (FilePath[i]!='/') && (FilePath[i]!='\\'))  i--;
      if (FilePath[i]=='.')  {
        FilePath[i+1] = char(0);
        strcat ( FilePath,newExt );
      } else  {
        strcat ( FilePath,"."    );
        strcat ( FilePath,newExt );
      }
      return FilePath;
    }

    bool FileExists ( cpstr FileName, PFile f )  {
    PFile g;
    bool  B;
      if (FileName)  {
        if (!f)  g = new File();
           else  g = f;
        g->assign ( FileName );
        B = g->exists();
        if (!f)  delete g;
        return B;
      } else
        return false;
    }


    //  ========================  File Class  ========================

    #define ARCH_NONE      0
    #define ARCH_GZIP      1
    #define ARCH_COMPRESS  2
    #define ARCH_ENFORCE   3

    File::File ( word BufSize )  {
      Buf_Size  = BufSize;
      BufLen    = 0;
      BufInc    = 1;
      EofFile   = false;
      hFile     = NULL;
      FName     = NULL;
      BufCnt    = 0;
      IOBuf     = NULL;
      IOSuccess = true;
      TextMode  = false;
      UniBin    = false;
      StdIO     = false;
      gzipIO    = ARCH_NONE;
      memIO     = false;
      ownBuf    = true;
    }

    File::~File()  {
      shut      ();
      FreeBuffer();
    }

    void  File::FreeBuffer ()  {
      if (IOBuf)  {
        if (ownBuf)  delete[] IOBuf;
        IOBuf = NULL;
      }
      if (FName)  {
        delete[] FName;
        FName = NULL;
      }
    }

    void  File::assign ( cpstr FileName, bool Text, bool UniB,
                         GZ_MODE gzMode )  {
    pstr p;

      shut();
      FreeBuffer();
      ownBuf = true;

      CreateCopy ( FName,FileName );
      StdIO = (!strcmp(FName,"stdin" )) ||
              (!strcmp(FName,"stdout")) ||
              (!strcmp(FName,"stderr"));

      if (StdIO)  TextMode = true;
            else  TextMode = Text;

      UniBin = UniB;

      gzipMode = gzMode;
      gzipIO   = ARCH_NONE;
      if ((gzipMode==GZM_ENFORCE) || (gzipMode==GZM_ENFORCE_GZIP))
        gzipIO = ARCH_GZIP;
      else if (gzipMode==GZM_ENFORCE_COMPRESS)
        gzipIO = ARCH_COMPRESS;
      else if (gzipMode==GZM_CHECK)  {
        p = LastOccurence ( FName,'.' );
        if (p)  {
          if (!strcmp(p,".gz"))      gzipIO = ARCH_GZIP;
          else if (!strcmp(p,".Z"))  gzipIO = ARCH_COMPRESS;
        }
      }

      memIO = false;

    }


    void  File::assign ( word poolSize, word sizeInc, pstr filePool )  {

      shut();

      IOBuf    = (pstr)filePool;
      BufLen   = poolSize;
      FLength  = poolSize;
      BufInc   = sizeInc;
      BufCnt   = 0;

      memIO    = true;
      ownBuf   = (IOBuf==NULL);
      gzipMode = GZM_NONE;
      gzipIO   = ARCH_NONE;

    }

    void  File::truncate ( long size )  {
    // call before reset/append
#ifndef _WIN32
      ::truncate ( FName,size );
#endif
    }

    void  File::takeFilePool ( pstr & filePool, word & fileSize )  {
      if (memIO)  {
        filePool = IOBuf;
        fileSize = FLength;
        IOBuf    = NULL;
        BufLen   = 0;
        BufCnt   = 0;
        FLength  = 0;
        ownBuf   = false;
      } else  {
        filePool = NULL;
        fileSize = 0;
      }
    }

    static pstr gzip_path       = pstr("gzip ");
    static pstr ungzip_path     = pstr("gzip -dc ");
    static pstr compress_path   = pstr("compress ");
    static pstr uncompress_path = pstr("uncompress -c ");

    void  SetGZIPPath ( pstr gzipPath, pstr ungzipPath )  {
      if (!gzipPath)  gzip_path = pstr("gzip ");
                else  gzip_path = gzipPath;
      if (!ungzipPath)  ungzip_path = pstr("gzip -d ");
                  else  ungzip_path = ungzipPath;
    }

    void  SetCompressPath ( pstr compressPath, pstr uncompressPath )  {
      if (!compressPath)  compress_path = pstr("compress ");
                    else  compress_path = compressPath;
      if (!uncompressPath)  uncompress_path = pstr("uncompress -c ");
                      else  uncompress_path = uncompressPath;
    }


    bool  File::reset ( bool ReadOnly, int retry )  {
    #ifndef _MSC_VER
    pstr p;
    int  i;
    #endif

      if (memIO)  {

        if (!IOBuf)  return false;
        BufCnt    = 0;
        IOSuccess = true;

      } else  {

        if (!FName)  return false;
        shut();
        BufLen = 0;
        BufCnt = 0;

        if (!strcmp(FName,"stdin"))  {

          hFile     = stdin;
          StdIO     = true;
          TextMode  = true;
          FLength   = 1;
          EofFile   = false;
          IOSuccess = true;

        } else  {

          StdIO = false;
          if (gzipIO==ARCH_GZIP)  {
    #ifndef _MSC_VER
            p = NULL;
            CreateConcat  ( p,ungzip_path,FName );
            for (i=0;(i<=retry) && (!hFile);i++)  {
              if (i>0)  sleep ( 1 );
              hFile = popen ( p,"r" );
            }
            if (p)  delete[] p;
    #endif

          } else if (gzipIO==ARCH_COMPRESS)  {
    #ifndef _MSC_VER
            p = NULL;
            CreateConcat  ( p,uncompress_path,FName );
            for (i=0;(i<=retry) && (!hFile);i++)  {
              if (i>0)  sleep ( 1 );
              hFile = popen ( p,"r" );
            }
            if (p)  delete[] p;
    #endif

          } else  {

    #ifndef _MSC_VER
            for (i=0;(i<=retry) && (!hFile);i++)  {
              if (i>0)  sleep ( 1 );
              if (TextMode)  {
                if (ReadOnly)  hFile = fopen ( FName,"rt"  );
                         else  hFile = fopen ( FName,"r+t" );
              } else  {
                if (ReadOnly)  hFile = fopen ( FName,"rb"  );
                         else  hFile = fopen ( FName,"r+b" );
              }
            }
    #endif

          }

          if (hFile)  {
            if (gzipIO==ARCH_NONE)  {
              fseek ( hFile,0L,SEEK_END );
              FLength = ftell ( hFile );
              fseek ( hFile,0L,SEEK_SET );
              EofFile = (FLength<=0);
            } else  {
              FLength = 1;
              EofFile = false;
            }
            IOSuccess = true;
          } else  {
            EofFile   = true;
            IOSuccess = false;
          }

        }

      }

      return  IOSuccess;

    }

    bool  File::rewrite()  {
    #ifndef _MSC_VER
    pstr p;
    #endif

      if (memIO)  {

        shut();
        if (IOBuf)  delete[] IOBuf;
        IOBuf = new char[BufLen];
        BufCnt    = 0;
        FLength   = 0;
        IOSuccess = true;
        ownBuf    = true;

      } else  {

        if (!FName)  return false;
        shut();
        BufLen = 0;
        BufCnt = 0;

        if (gzipIO==ARCH_GZIP)  {
    #ifndef _MSC_VER
          p = NULL;
          CreateConcat  ( p,gzip_path,pstr(" > "),FName );
          hFile = popen ( p,"w" );
          if (p)  delete[] p;
    #else
          hFile = NULL;
    #endif
          StdIO = false;
        } else if (gzipIO==ARCH_COMPRESS)  {
    #ifndef _MSC_VER
          p = NULL;
          CreateConcat  ( p,compress_path,pstr(" > "),FName );
          hFile = popen ( p,"w" );
          if (p)  delete[] p;
    #else
          hFile = NULL;
    #endif
          StdIO = false;
        } else if (!TextMode)  {
          hFile = fopen ( FName,"w+b" );
          StdIO = false;
        } else if (!strcmp(FName,"stdout"))  {
          hFile = stdout;
          StdIO = true;
        } else if (!strcmp(FName,"stderr")) {
          hFile = stderr;
          StdIO = true;
        } else  {
          hFile = fopen ( FName,"w+t" );
          StdIO = false;
        }

        FLength   = 0;
        IOSuccess = (hFile!=NULL);

      }

      return  IOSuccess;

    }


    bool  File::append()  {
    #ifndef _MSC_VER
    pstr p;
    #endif

      if (memIO)  {

        if (!IOBuf)  {
          IOBuf  = new char[BufLen];
          BufCnt = 0;
          ownBuf = true;
        }
        FLength   = BufCnt;
        IOSuccess = true;

      } else  {

        if (!FName)  return false;

        shut();
        BufLen  = 0;
        BufCnt  = 0;
        if (gzipIO==ARCH_GZIP)  {
    #ifndef _MSC_VER
          p = NULL;
          CreateConcat  ( p,gzip_path,pstr(" >> "),FName );
          hFile = popen ( p,"w" );
          if (p)  delete[] p;
    #else
          hFile = NULL;
    #endif
          StdIO = false;
        } else if (gzipIO==ARCH_COMPRESS)  {
    #ifndef _MSC_VER
          p = NULL;
          CreateConcat  ( p,compress_path,pstr(" >> "),FName );
          hFile = popen ( p,"w" );
          if (p)  delete[] p;
    #else
          hFile = NULL;
    #endif
          StdIO = false;
        } else if (!TextMode)  {
          hFile = fopen ( FName,"ab" );
          StdIO = false;
        } else if (!strcmp(FName,"stdout"))  {
          hFile = stdout;
          StdIO = true;
        } else if (!strcmp(FName,"stderr")) {
          hFile = stderr;
          StdIO = true;
        } else  {
          hFile = fopen ( FName,"at" );
          StdIO = false;
        }

        FLength = 0;
        IOSuccess = hFile!=NULL;

      }

      return  IOSuccess;

    }

    bool  File::erase()  {
      if (!FName)  return false;
      shut();
      if (!StdIO)  {
        BufLen = 0;
        BufCnt = 0;
        if (FName)
          IOSuccess = (remove(FName)==0);
        FLength = 0;
      } else
        IOSuccess = true;
      return  IOSuccess;
    }

    bool  File::exists()  {

      if (memIO)  {

        IOSuccess = (IOBuf!=NULL);

      } else  {

        if (!FName)  return false;
        shut();
        if (!StdIO)  {
          hFile     = fopen ( FName,"r" );
          IOSuccess = (hFile!=NULL);
          BufLen    = 0;
          BufCnt    = 0;
          FLength   = 0;
          if (hFile)  fclose ( hFile );
        } else
          IOSuccess = true;
        hFile = NULL;
      }

      return  IOSuccess;

    }

    bool  File::parse ( cpstr FileName )  {
    UNUSED_ARGUMENT(FileName);
      return true;
    }

    bool  File::rename ( cpstr NewFileName )  {
      if (!FName)  return false;
      shut();
      if (!StdIO)
        IOSuccess = (::rename(FName,NewFileName)==0);
      if (IOSuccess)  assign ( NewFileName,TextMode,UniBin,gzipMode );
      return  IOSuccess;
    }

    long  File::Position()  {
    // do not use on text files

      if (memIO)  return BufCnt;

      if (hFile==NULL)  return 0L;
      return  ftell ( hFile );

    }

    bool  File::seek ( long Position )  {
    // do not use on text files
      if (memIO)  {
        if (Position<=(long)BufLen)  {
          BufCnt    = Position;
          IOSuccess = true;
        } else
          IOSuccess = false;
        return IOSuccess;
      } else if (hFile==NULL)
        return false;
      else if (!StdIO)  {
        IOSuccess = fseek(hFile,Position,SEEK_SET)==0;
        return IOSuccess;
      } else
        return true;
    }

    bool  File::FileEnd()  {

      if (memIO)  return ((long)BufCnt>=FLength);

      if (TextMode)  {
        if (EofFile || ((!hFile) && (!StdIO)))
           return true;
        if (feof(hFile)==0)
           return false;
        return true;
      }

      return  EofFile && (BufLen==0);

    }

    void  File::flush ()  {

      if (hFile!=NULL)  {
        if (!StdIO)  {
    #ifndef _MSC_VER
          if (gzipIO==ARCH_NONE)
            fflush ( hFile );
    #else
          fflush ( hFile );
    #endif
        }
      }
    }

    void  File::shut ()  {
    
      if (memIO)
        FreeBuffer();

      else if (hFile!=NULL)  {
        if (!StdIO)  {
    #ifndef _MSC_VER
          if (gzipIO!=ARCH_NONE)  pclose ( hFile );
                            else  fclose ( hFile );
    #else
          fclose ( hFile );
    #endif
        }
        hFile = NULL;
      }

    }

    bool File::isOpen()  {
      if (memIO)  return (IOBuf!=NULL);
      return (hFile!=NULL);
    }

    word  File::ReadLine ( pstr Line, word MaxLen )  {
    word     LCnt;
    int      Done;
    bool  HSuccess = IOSuccess;

      if (memIO)  {

        LCnt = 0;
        while (((long)BufCnt<FLength) &&
               (LCnt<MaxLen-1)        &&
               (IOBuf[BufCnt]!='\r')  &&
               (IOBuf[BufCnt]!='\n') )
          Line[LCnt++] = IOBuf[BufCnt++];
        Line[LCnt] = char(0);

        while (((long)BufCnt<FLength) &&
               (IOBuf[BufCnt]!='\r')  &&
               (IOBuf[BufCnt]!='\n') )
          BufCnt++;

        if ((long)BufCnt<FLength)  {
          if (IOBuf[BufCnt]=='\r')  {
            BufCnt++;
            if ((long)BufCnt<FLength)  {
              if (IOBuf[BufCnt]=='\n')  BufCnt++;
            }
          } else if (IOBuf[BufCnt]=='\n')  {
            BufCnt++;
            if ((long)BufCnt<FLength)  {
              if (IOBuf[BufCnt]=='\r')  BufCnt++;
            }
          }
        }

        return LCnt;

      } else  {

        if ((!hFile) && (!StdIO))  {
          Line[0]   = char(0);
          EofFile   = true;
          BufLen    = 0;
          IOSuccess = false;
          return 0;
        }
        if (TextMode)  {
          Line[0] = char(0);
          if (fgets(Line,MaxLen,hFile))  {
            LCnt = strlen(Line);
            while (LCnt>0)  {
              if ((Line[LCnt-1]!='\n') &&
                  (Line[LCnt-1]!='\r'))  break;
              Line[LCnt-1] = char(0);
              LCnt--;
            }
          } else
            LCnt = 0;
          return LCnt;
        } else  {
          if (IOBuf==NULL)  {
            IOBuf     = new char[Buf_Size];
            BufLen    = ReadFile ( IOBuf,Buf_Size );
            IOSuccess = HSuccess;
            BufCnt    = 0;
            ownBuf    = true;
          }
          LCnt = 0;
          do {
            while ((BufCnt<BufLen)       &&
                   (LCnt<MaxLen-1)       &&
                   (IOBuf[BufCnt]!='\r') &&
                   (IOBuf[BufCnt]!='\n') )
              Line[LCnt++] = IOBuf[BufCnt++];
            if (BufCnt>=BufLen)  {
              HSuccess  = IOSuccess;
              BufLen    = ReadFile ( IOBuf,Buf_Size );
              IOSuccess = HSuccess;
              BufCnt    = 0;
            }
            if (IOBuf[BufCnt]=='\r')      Done = 1;
            else if (IOBuf[BufCnt]=='\n') Done = 2;
                                     else Done = 0;
            if (Done)  BufCnt++;
            if (BufCnt>=BufLen)  {
              HSuccess  = IOSuccess;
              BufLen    = ReadFile ( IOBuf,Buf_Size );
              IOSuccess = HSuccess;
              BufCnt    = 0;
            }
            if (BufLen>0)  {
              if (((Done==2) && (IOBuf[BufCnt]=='\r')) ||
                  ((Done==1) && (IOBuf[BufCnt]=='\n')))
                BufCnt++;
            }
          } while ((!Done) && (LCnt<MaxLen-1) && (BufLen>0));
          Line[LCnt] = char(0);
          return LCnt;
        }

      }

    }

    word  File::ReadNonBlankLine ( pstr S, word MaxLen )  {
    word  i,j;
      do  {
        j = ReadLine ( S,MaxLen );
        i = 0;
        while ((i<j) && (S[i]==' ')) i++;
      } while ((i>=j) && (!FileEnd()));
      if (i>=j)  {
        S[0] = char(0);
        j    = 0;
      }
      return  j;
    }


    bool  File::WriteLine ( cpstr Line )  {
      if ((!memIO) && TextMode)  {
        if (hFile==NULL)  return false;
        fputs ( Line,hFile );
    //    return (fputs(NEWLINE,hFile)>=0);
        return (fputs("\n",hFile)>=0);
      } else  {
        if (WriteFile(Line,strlen(Line)))
              return  WriteFile ( (void *)NEWLINE,strlen(NEWLINE) );
        else  return  false;
      }
    }

    bool  File::Write ( cpstr Line )  {
      if ((!memIO) && TextMode)  {
        if (hFile==NULL)  return false;
        return (fputs(Line,hFile)>=0);
      } else
        return WriteFile(Line,strlen(Line));
    }

    bool  File::Write ( realtype V, int length )  {
    char N[50];
      sprintf ( N,"%-.*g",length,V );
      if ((!memIO) && TextMode)  {
        if (hFile==NULL)  return false;
        return (fputs(N,hFile)>=0);
      } else
        return WriteFile(N,strlen(N));
    }

    bool  File::Write ( int iV, int length )  {
    char N[50];
      sprintf ( N,"%*i",length,iV );
      if ((!memIO) && TextMode)  {
        if (hFile==NULL)  return false;
        return (fputs(N,hFile)>=0);
      } else
        return WriteFile(N,strlen(N));
    }

    bool  File::LF()  {
      if ((!memIO) && TextMode)  {
        if (hFile==NULL)  return false;
    //    return (fputs(NEWLINE,hFile)>=0);
        return (fputs("\n",hFile)>=0);
      } else
        return WriteFile ( (void *)NEWLINE,strlen(NEWLINE) );
    }

    bool  File::WriteDataLine ( realtype X, realtype Y, int length )  {
      Write ( pstr("   ") );
      Write ( X,length    );
      Write ( pstr("   ") );
      Write ( Y,length    );
      return LF();
    }

    bool  File::WriteParameter ( cpstr S, realtype X,
                                 int ParColumn, int length )  {
    int  l=strlen(S);
      if ((!memIO) && TextMode)  {
        fputs ( S,hFile );
        while (l<ParColumn)  {
          fputs ( " ",hFile );
          l++;
        }
      } else  {
        WriteFile ( S,l );
        while (l<ParColumn)  {
          WriteFile ( (void *)" ",1 );
          l++;
        }
      }
      Write ( X,length );
      return LF();
    }

    bool  File::WriteParameters ( cpstr S, int n_X, rvector X,
                                  int ParColumn, int length )  {
    int  i;
    int  l=strlen(S);
      if ((!memIO) && TextMode)  {
        fputs ( S,hFile );
        while (l<ParColumn)  {
          fputs ( " ",hFile );
          l++;
        }
      } else  {
        WriteFile ( S,l );
        while (l<ParColumn)  {
          WriteFile ( (void *)" ",1 );
          l++;
        }
      }
      for (i=0;i<n_X;i++)  {
        Write ( X[i],length );
        if (i!=n_X-1)  WriteFile ( (void *)", ",2 );
      }
      return LF();
    }

    bool  File::ReadParameter  ( pstr S, realtype & X,
                                 int ParColumn )  {
      ReadLine ( S );
      if ((int)strlen(S)>ParColumn)  {
    //    X = atof ( &(S[ParColumn]) );
        X = GetNumber ( &(S[ParColumn]) );
        return true;
      } else  {
        X = 0.0;
        return false;
      }
    }

    bool  File::ReadParameters ( pstr S, int & n_X, rvector X,
                                 int MaxLen, int ParColumn )  {
    pstr S1,S2;
      ReadLine ( S,MaxLen );
      if ((int)strlen(S)>ParColumn)  {
        n_X = 0;
        S2 = &(S[ParColumn]);
        S1 = S2;
        while (*S1!=char(0))  {
          if (*S1==',')  *S1 = ' ';
          S1++;
        }
        while (*S2!=char(0))  {
          S1 = S2;
          X[n_X] = strtod ( S1,&S2 );
          n_X++;
          while ((*S2!=char(0)) && (*S2==' '))  S2++;
        }
        return true;
      } else  {
        n_X  = 0;
        X[0] = 0.0;
        return false;
      }
    }

    bool  File::ReadParameter  ( pstr S, int & X,
                                 int ParColumn )  {
    realtype  V;
      if (ReadParameter(S,V,ParColumn))  {
        X = mround(V);
        return true;
      } else  {
        X = 0;
        return false;
      }
    }

    bool  File::CreateWrite ( cpstr Line )  {
    wordUniBin wUB;
    word       i;
      if (UniBin)  {
        if (Line)  {
          i = strlen(Line)+1;
          word2UniBin ( i,wUB );
          if (WriteFile(wUB,sizeof(wordUniBin)))
                 return WriteFile ( Line,i );
          else   return false;
        } else  {
          i = 0;
          word2UniBin ( i,wUB );
          return WriteFile ( wUB,sizeof(wordUniBin) );
        }
      } else  {
        if (Line)  {
          i = strlen(Line)+1;
          if (WriteFile(&i,sizeof(i)))
                 return WriteFile ( Line,i );
          else   return false;
        } else  {
          i = 0;
          return WriteFile ( &i,sizeof(i) );
        }
      }
    }

    #define _max_dyn_string_len  1073741824

    word  File::CreateRead ( pstr & Line )  {
    wordUniBin wUB;
    word       i;
    //unsigned short int i;
      if (Line)  {
        delete[] Line;
        Line = NULL;
      }
      if (UniBin)  {
        ReadFile    ( wUB,sizeof(wordUniBin) );
        UniBin2word ( wUB,i );
      } else
        ReadFile ( &i,sizeof(i) );
      if ((i>0) && (i<_max_dyn_string_len))  {
        Line = new char[i];
        ReadFile ( Line,i );
      }
      return i;
    }

    bool  File::WriteTerLine ( cpstr Line, bool longLine )  {
    wordUniBin wUB;
    word       ll;
    byte       sl;
    bool    B;
      if (Line)  ll = strlen(Line);
           else  ll = 0;
      if (!longLine)  {
        sl = byte(ll);
        B  = WriteFile ( &sl,sizeof(sl) );
      } else if (UniBin)  {
        word2UniBin ( ll,wUB );
        B = WriteFile ( wUB,sizeof(wordUniBin) );
      } else
        B = WriteFile ( &ll,sizeof(ll) );
      if (B && (ll>0))  B = WriteFile ( Line,ll );
      return B;
    }

    word  File::ReadTerLine ( pstr Line, bool longLine )  {
    wordUniBin wUB;
    word       ll;
    byte       sl;
      if (!longLine)  {
        ReadFile ( &sl,sizeof(sl) );
        ll = sl;
      } else if (UniBin)  {
        ReadFile    ( wUB,sizeof(wordUniBin) );
        UniBin2word ( wUB,ll );
      } else
        ReadFile ( &ll,sizeof(ll) );
      if (ll>0)  ReadFile ( Line,ll );
      Line[ll] = char(0);
      return ll+1;
    }

    word  File::ReadFile ( void * Buffer, word Count )  {
    word  Cnt;
      if (memIO)  {
        Cnt       = WMin(Count,FLength-BufCnt);
        if (Cnt>0)  {
          memcpy ( Buffer,&(IOBuf[BufCnt]),Cnt );
          BufCnt += Cnt;
        }
        IOSuccess = (Cnt==Count);
        EofFile   = ((Cnt<Count) || ((long)BufCnt>=FLength));
        return  Cnt;
      } else if (hFile)  {
        Cnt       = (word)fread ( Buffer,1,Count,hFile );
        EofFile   = (Cnt<Count) ||
                     ((gzipIO==ARCH_NONE) && (Position()==FLength));
        IOSuccess = (Cnt==Count);
        return  Cnt;
      } else
        return  0;
    }

    bool  File::WriteFile ( const void * Buffer, word Count )  {
    pstr IOB;
    word Cnt;
    long Pos;

      if (memIO)  {

        Cnt = BufCnt + Count;
        if (Cnt>BufLen)  {
          Cnt += BufInc;
          IOB = new char[Cnt];
          if (IOBuf)  {
            memcpy ( IOB,IOBuf,BufCnt );
            delete[] IOBuf;
          }
          IOBuf  = IOB;
          BufLen = Cnt;
          ownBuf = true;
        }
        memcpy ( &(IOBuf[BufCnt]),Buffer,Count );
        BufCnt += Count;
        FLength = BufCnt;
        IOSuccess = true;

      } else  {

        if (hFile==NULL)  return false;
        Cnt = (word)fwrite ( Buffer,1,Count,hFile );
        Pos = Position();
        if (Pos>FLength)  FLength = Pos;
        IOSuccess = Cnt==Count;

      }

      return IOSuccess;

    }

    bool  File::WriteReal ( realtype * V )  {
    realUniBin  rUB;
      if (UniBin)  {
        real2UniBin ( *V,rUB );
        return WriteFile ( rUB,sizeof(realUniBin) );
      } else
        return WriteFile ( V,sizeof(realtype) );
    }

    bool  File::WriteFloat ( realtype * V )  {
    floatUniBin fUB;
    float       fV;
      if (UniBin)  {
        float2UniBin ( *V,fUB );
        return WriteFile ( fUB,sizeof(floatUniBin) );
      } else  {
        fV = (float)*V;
        return WriteFile ( &fV,sizeof(float) );
      }
    }

    bool  File::WriteInt ( int * I )  {
    intUniBin  iUB;
      if (UniBin)  {
        int2UniBin ( *I,iUB );
        return WriteFile ( iUB,sizeof(intUniBin) );
      } else
        return WriteFile ( I,sizeof(int) );
    }

    bool  File::WriteShort ( short * S )  {
    shortUniBin  sUB;
      if (UniBin)  {
        short2UniBin ( *S,sUB );
        return WriteFile ( sUB,sizeof(shortUniBin) );
      } else
        return WriteFile ( S,sizeof(short) );
    }

    bool  File::WriteLong ( long * L )  {
    longUniBin  lUB;
      if (UniBin)  {
        long2UniBin ( *L,lUB );
        return WriteFile ( lUB,sizeof(longUniBin) );
      } else
        return WriteFile ( L,sizeof(long) );
    }

    bool  File::WriteBool ( bool * B )  {
    intUniBin  iUB;
    int        k;
      if (UniBin)  {
        if (*B)  k = 1;
           else  k = 0;
        int2UniBin ( k,iUB );
        return WriteFile ( iUB,sizeof(intUniBin) );
      } else
        return WriteFile ( B,sizeof(bool) );
    }

    bool  File::WriteByte ( byte * B )  {
      return WriteFile ( B,sizeof(byte) );
    }

    bool  File::WriteWord ( word * W )  {
    wordUniBin  wUB;
      if (UniBin)  {
        word2UniBin ( *W,wUB );
        return WriteFile ( wUB,sizeof(wordUniBin) );
      } else
        return WriteFile ( W,sizeof(word) );
    }

    bool  File::ReadReal ( realtype * V )  {
    realUniBin  rUB;
      if (UniBin)  {
        if (ReadFile(rUB,sizeof(realUniBin))==sizeof(realUniBin))  {
          UniBin2real ( rUB,*V );
          return true;
        } else
          return false;
      } else
        return ( ReadFile(V,sizeof(realtype))==sizeof(realtype) );
    }

    bool  File::ReadFloat ( realtype * V )  {
    floatUniBin  fUB;
    float        fV;
      if (UniBin)  {
        if (ReadFile(fUB,sizeof(floatUniBin))==sizeof(floatUniBin))  {
          UniBin2float ( fUB,*V );
          return true;
        }
      } else if (ReadFile(&fV,sizeof(float))==sizeof(float))  {
        *V = fV;
        return true;
      }
      return false;
    }

    bool  File::ReadInt ( int * I )  {
    intUniBin  iUB;
      if (UniBin)  {
        if (ReadFile(iUB,sizeof(intUniBin))==sizeof(intUniBin))  {
          UniBin2int ( iUB,*I );
          return true;
        } else
          return false;
      } else
        return ( ReadFile(I,sizeof(int))==sizeof(int) );
    }

    bool  File::ReadShort ( short * S )  {
    shortUniBin  sUB;
      if (UniBin)  {
        if (ReadFile(sUB,sizeof(shortUniBin))==sizeof(shortUniBin))  {
          UniBin2short ( sUB,*S );
          return true;
        } else
          return false;
      } else
        return ( ReadFile(S,sizeof(short))==sizeof(short) );
    }

    bool  File::ReadLong ( long * L )  {
    longUniBin  lUB;
      if (UniBin)  {
        if (ReadFile(lUB,sizeof(longUniBin))==sizeof(longUniBin))  {
          UniBin2long ( lUB,*L );
          return true;
        } else
          return false;
      } else
        return ( ReadFile(L,sizeof(long))==sizeof(long) );
    }

    bool  File::ReadBool ( bool * B )  {
    intUniBin  iUB;
    int        k;
      if (UniBin)  {
        if (ReadFile(iUB,sizeof(intUniBin))==sizeof(intUniBin))  {
          UniBin2int ( iUB,k );
          *B = (k!=0);
          return true;
        } else
          return false;
      } else
        return ( ReadFile(B,sizeof(bool))==sizeof(bool) );
    }

    bool  File::ReadByte ( byte * B )  {
      return ( ReadFile(B,sizeof(byte))==sizeof(byte) );
    }

    bool  File::ReadWord ( word * W )  {
    wordUniBin  wUB;
      if (UniBin)  {
        if (ReadFile(wUB,sizeof(wordUniBin))==sizeof(wordUniBin))  {
          UniBin2word ( wUB,*W );
          return true;
        } else
          return false;
      } else
        return ( ReadFile(W,sizeof(word))==sizeof(word) );
    }

    bool  File::AddReal ( realtype * V )  {
    realtype x;
      if (ReadReal(&x))  {
        *V += x;
        return true;
      }
      return false;
    }

    bool  File::AddFloat ( realtype * V )  {
    realtype x;
      if (ReadFloat(&x))  {
        *V += x;
        return true;
      }
      return false;
    }


    bool  File::AddInt ( int * I )  {
    int k;
      if (ReadInt(&k))  {
        *I += k;
        return true;
      }
      return false;
    }

    bool  File::AddShort ( short * S )  {
    short k;
      if (ReadShort(&k))  {
        *S += k;
        return true;
      }
      return false;
    }

    bool  File::AddLong ( long * L )  {
    long k;
      if (ReadLong(&k))  {
        *L += k;
        return true;
      }
      return false;
    }

    bool  File::AddByte ( byte * B )  {
    byte k;
      if (ReadByte(&k))  {
        *B += k;
        return true;
      }
      return false;
    }

    bool  File::AddWord ( word * W )  {
    word k;
      if (ReadWord(&k))  {
        *W += k;
        return true;
      }
      return false;
    }

    bool  File::WriteVector ( rvector V, int len, int Shift )  {
    intUniBin  iUB;
    realUniBin rUB;
    int        i;
    int        l = len;
      if (V==NULL)  l = 0;
      if (UniBin)  {
        int2UniBin ( l,iUB );
        WriteFile ( iUB,sizeof(intUniBin) );
        for (i=0;i<len;i++)  {
          real2UniBin ( V[Shift+i],rUB );
          WriteFile   ( rUB,sizeof(realUniBin) );
        }
      } else  {
        WriteFile ( &l,sizeof(l) );
        if (l>0)  WriteFile ( &(V[Shift]),sizeof(realtype)*l );
      }
      return IOSuccess;
    }

    bool  File::WriteVector ( ivector iV, int len, int Shift )  {
    intUniBin  iUB;
    int        i;
    int        l = len;
      if (iV==NULL)  l = 0;
      if (UniBin)  {
        int2UniBin ( l,iUB );
        WriteFile ( iUB,sizeof(intUniBin) );
        for (i=0;i<len;i++)  {
          int2UniBin ( iV[Shift+i],iUB );
          WriteFile  ( iUB,sizeof(intUniBin) );
        }
      } else  {
        WriteFile ( &l,sizeof(l) );
        if (l>0)  WriteFile ( &(iV[Shift]),sizeof(int)*l );
      }
      return IOSuccess;
    }

    bool  File::WriteVector ( lvector lV, int len, int Shift )  {
    intUniBin  iUB;
    longUniBin lUB;
    int        i;
    int        l = len;
      if (lV==NULL)  l = 0;
      if (UniBin)  {
        int2UniBin ( l,iUB );
        WriteFile ( iUB,sizeof(intUniBin) );
        for (i=0;i<len;i++)  {
          long2UniBin ( lV[Shift+i],lUB );
          WriteFile   ( lUB,sizeof(longUniBin) );
        }
      } else  {
        WriteFile ( &l,sizeof(l) );
        if (l>0)  WriteFile ( &(lV[Shift]),sizeof(long)*l );
      }
      return IOSuccess;
    }

    bool  File::WriteVector ( bvector B, int len, int Shift )  {
    intUniBin iUB;
    int       l = len;
      if (B==NULL)  l = 0;
      if (UniBin)  {
        int2UniBin ( l,iUB );
        WriteFile  ( iUB,sizeof(intUniBin) );
      } else
        WriteFile ( &l,sizeof(l) );
      if (l>0)  WriteFile ( &(B[Shift]),sizeof(byte)*l );
      return IOSuccess;
    }

    bool  File::ReadVector ( rvector V, int maxlen, int Shift )  {
    intUniBin  iUB;
    realUniBin rUB;
    int        i,l,ll;
    realtype   B;
      if (UniBin)  {
        ReadFile   ( iUB,sizeof(intUniBin) );
        UniBin2int ( iUB,l );
        if (IOSuccess && (l>0))  {
          ll = IMin(l,maxlen);
          if (V)
            for (i=0;i<=ll;i++)  {
              ReadFile    ( rUB,sizeof(realUniBin) );
              UniBin2real ( rUB,V[Shift+i] );
            }
          for (i=ll+1;i<=l;i++)
            ReadFile ( rUB,sizeof(realUniBin) );
        }
      } else  {
        ReadFile ( &l,sizeof(l) );
        if (IOSuccess && (l>0))  {
          ll = IMin(l,maxlen);
          if (V)  ReadFile ( &(V[Shift]),sizeof(realtype)*ll );
          for (i=ll+1;i<=l;i++)  ReadFile ( &B,sizeof(B) );
        }
      }
      return IOSuccess;
    }

    bool  File::ReadVector ( ivector iV, int maxlen, int Shift )  {
    intUniBin  iUB;
    int        i,l,ll,iB;
      if (UniBin)  {
        ReadFile   ( iUB,sizeof(intUniBin) );
        UniBin2int ( iUB,l );
        if (IOSuccess && (l>0))  {
          ll = IMin(l,maxlen);
          if (iV)
            for (i=0;i<=ll;i++)  {
              ReadFile   ( iUB,sizeof(intUniBin) );
              UniBin2int ( iUB,iV[Shift+i]       );
            }
          for (i=ll+1;i<=l;i++)
            ReadFile ( iUB,sizeof(intUniBin) );
        }
      } else  {
        ReadFile ( &l,sizeof(l) );
        if (IOSuccess && (l>0))  {
          ll = IMin(l,maxlen);
          if (iV)  ReadFile ( &(iV[Shift]),sizeof(int)*ll );
          for (i=ll+1;i<=l;i++)  ReadFile ( &iB,sizeof(iB) );
        }
      }
      return IOSuccess;
    }

    bool  File::ReadVector ( lvector lV, int maxlen, int Shift )  {
    intUniBin  iUB;
    longUniBin lUB;
    int        i,l,ll;
    long       lB;
      if (UniBin)  {
        ReadFile   ( iUB,sizeof(intUniBin) );
        UniBin2int ( iUB,l );
        if (IOSuccess && (l>0))  {
          ll = IMin(l,maxlen);
          if (lV)
            for (i=0;i<=ll;i++)  {
              ReadFile    ( lUB,sizeof(longUniBin) );
              UniBin2long ( lUB,lV[Shift+i]       );
            }
          for (i=ll+1;i<=l;i++)
            ReadFile ( lUB,sizeof(longUniBin) );
        }
      } else  {
        ReadFile ( &l,sizeof(l) );
        if (IOSuccess && (l>0))  {
          ll = IMin(l,maxlen);
          if (lV)  ReadFile ( &(lV[Shift]),sizeof(long)*ll );
          for (i=ll+1;i<=l;i++)  ReadFile ( &lB,sizeof(lB) );
        }
      }
      return IOSuccess;
    }

    bool  File::ReadVector ( bvector B, int maxlen, int Shift )  {
    intUniBin  iUB;
    int        i,l,ll;
    byte       t;
      if (UniBin)  {
        ReadFile   ( iUB,sizeof(intUniBin) );
        UniBin2int ( iUB,l );
      } else
        ReadFile ( &l,sizeof(l) );
      if (IOSuccess && (l>0))  {
        ll = IMin(l,maxlen);
        if (B)  ReadFile ( &(B[Shift]),sizeof(byte)*ll );
        for (i=ll+1;i<=l;i++)  ReadFile ( &t,sizeof(t) );
      }
      return IOSuccess;
    }

    bool  File::CreateReadVector ( rvector & V, int & len,
                                       int Shift )  {
    intUniBin  iUB;
    realUniBin rUB;
    int        i;
    realtype   B;
      FreeVectorMemory ( V,Shift );
      if (UniBin)  {
        ReadFile   ( iUB,sizeof(intUniBin) );
        UniBin2int ( iUB,len );
        if (IOSuccess && (len>0))  {
          GetVectorMemory ( V,len,Shift );
          if (V)
            for (i=0;i<len;i++)  {
              ReadFile    ( rUB,sizeof(realUniBin) );
              UniBin2real ( rUB,V[Shift+i] );
            }
          else for (i=0;i<len;i++)
            ReadFile ( rUB,sizeof(realUniBin) );
        }
      } else  {
        ReadFile ( &len,sizeof(len) );
        if (IOSuccess && (len>0))  {
          GetVectorMemory ( V,len,Shift );
          if (V)  ReadFile ( &(V[Shift]),sizeof(realtype)*len );
            else for (i=0;i<len;i++)
                   ReadFile ( &B,sizeof(B) );
        }
      }
      return IOSuccess;
    }

    bool  File::CreateReadVector ( rvector & V, int Shift )  {
    int len;
      return CreateReadVector ( V,len,Shift );
    }

    bool  File::CreateReadVector ( ivector & iV, int & len,
                                       int Shift )  {
    intUniBin  iUB;
    int        i,iB;
      FreeVectorMemory ( iV,Shift );
      if (UniBin)  {
        ReadFile   ( iUB,sizeof(intUniBin) );
        UniBin2int ( iUB,len );
        if (IOSuccess && (len>0))  {
          GetVectorMemory ( iV,len,Shift );
          if (iV)
            for (i=0;i<len;i++)  {
              ReadFile   ( iUB,sizeof(intUniBin) );
              UniBin2int ( iUB,iV[Shift+i]       );
            }
          else for (i=0;i<len;i++)
            ReadFile ( iUB,sizeof(intUniBin) );
        }
      } else  {
        ReadFile ( &len,sizeof(len) );
        if (IOSuccess && (len>0))  {
          GetVectorMemory ( iV,len,Shift );
          if (iV)  ReadFile ( &(iV[Shift]),sizeof(int)*len );
             else for (i=0;i<len;i++)
                    ReadFile ( &iB,sizeof(iB) );
        }
      }
      return IOSuccess;
    }

    bool  File::CreateReadVector ( ivector & iV, int Shift )  {
    int len;
      return CreateReadVector ( iV,len,Shift );
    }

    bool  File::CreateReadVector ( lvector & lV, int & len,
                                       int Shift )  {
    intUniBin  iUB;
    longUniBin lUB;
    int        i;
    long       lB;
      FreeVectorMemory ( lV,Shift );
      if (UniBin)  {
        ReadFile   ( iUB,sizeof(intUniBin) );
        UniBin2int ( iUB,len );
        if (IOSuccess && (len>0))  {
          GetVectorMemory ( lV,len,Shift );
          if (lV)
            for (i=0;i<len;i++)  {
              ReadFile    ( lUB,sizeof(intUniBin) );
              UniBin2long ( lUB,lV[Shift+i]       );
            }
          else for (i=0;i<len;i++)
            ReadFile ( lUB,sizeof(longUniBin) );
        }
      } else  {
        ReadFile ( &len,sizeof(len) );
        if (IOSuccess && (len>0))  {
          GetVectorMemory ( lV,len,Shift );
          if (lV) ReadFile ( &(lV[Shift]),sizeof(long)*len );
             else for (i=0;i<len;i++)
                    ReadFile ( &lB,sizeof(lB) );
        }
      }
      return IOSuccess;
    }

    bool  File::CreateReadVector ( lvector & lV, int Shift )  {
    int len;
      return CreateReadVector ( lV,len,Shift );
    }


    bool  File::CreateReadVector ( bvector & B, int & len,
                                       int Shift )  {
    intUniBin  iUB;
    int        i;
    byte       t;
      FreeVectorMemory ( B,Shift );
      if (UniBin)  {
        ReadFile   ( iUB,sizeof(intUniBin) );
        UniBin2int ( iUB,len );
      } else
        ReadFile ( &len,sizeof(len) );
      if (IOSuccess && (len>0))  {
        GetVectorMemory ( B,len,Shift );
        if (B)  ReadFile ( &(B[Shift]),sizeof(byte)*len );
          else for (i=0;i<len;i++)
                 ReadFile ( &t,sizeof(t) );
      }
      return IOSuccess;
    }

    bool  File::CreateReadVector ( bvector & B, int Shift )  {
    int len;
      return CreateReadVector ( B,len,Shift );
    }


    bool  File::WriteMatrix ( rmatrix & A, int N, int M,
                                  int ShiftN, int ShiftM )  {
    intUniBin  iUB;
    realUniBin rUB;
    int        i,j;
      if (UniBin)  {
        if (!A)  {
          i = 0;
          int2UniBin ( i,iUB );
          WriteFile  ( iUB,sizeof(intUniBin) );
        } else  {
          int2UniBin ( N,iUB );
          WriteFile  ( iUB,sizeof(intUniBin) );
          int2UniBin ( M,iUB );
          WriteFile  ( iUB,sizeof(intUniBin) );
          for (i=0;i<N;i++)
            for (j=0;j<M;j++)  {
              real2UniBin ( A[ShiftN+i][ShiftM+j],rUB );
              WriteFile   ( rUB,sizeof(realUniBin) );
            }
        }
      } else if (!A) {
        i = 0;
        WriteFile ( &i,sizeof(i) );
      } else  {
        WriteFile ( &N,sizeof(N) );
        WriteFile ( &M,sizeof(M) );
        for (i=0;i<N;i++)
          WriteFile ( &(A[ShiftN][ShiftM]),sizeof(realtype)*M );
      }
      return IOSuccess;
    }

    bool  File::CreateReadMatrix ( rmatrix & A, int ShiftN,
                                       int ShiftM )  {
    int N,M;
      return CreateReadMatrix ( A,N,M,ShiftN,ShiftM );
    }

    bool  File::CreateReadMatrix ( rmatrix & A, int & N, int & M,
                                       int ShiftN, int ShiftM )  {
    intUniBin  iUB;
    realUniBin rUB;
    int        i,j;
      FreeMatrixMemory ( A,N,ShiftN,ShiftM );
      if (UniBin)  {
        ReadFile   ( iUB,sizeof(intUniBin) );
        UniBin2int ( iUB,N );
        if (IOSuccess && (N>0))  {
          ReadFile   ( iUB,sizeof(intUniBin) );
          UniBin2int ( iUB,M );
          if (IOSuccess && (M>0))  {
            GetMatrixMemory ( A,N,M,ShiftN,ShiftM );
            for (i=0;i<N;i++)
              for (j=0;j<M;j++)  {
                ReadFile    ( rUB,sizeof(realUniBin) );
                UniBin2real ( rUB,A[ShiftN+i][ShiftM+j] );
              }
          }
        }
      } else  {
        ReadFile ( &N,sizeof(N) );
        if (N>0)  {
          ReadFile ( &M,sizeof(M) );
          if (M>0)  {
            GetMatrixMemory ( A,N,M,ShiftN,ShiftM );
            for (i=0;i<N;i++)
              ReadFile ( &(A[ShiftN][ShiftM]),sizeof(realtype)*M );
          }
        }
      }
      return IOSuccess;
    }


    bool  File::WriteColumns ( rvector X, rvector Y, rvector Z,
                                   int len, int Shift, int MLength )  {
    //   WriteColumns writes data stored in X, Y and Z in the form
    // of columns, adding a blank line in the end. If Z (or Z and Y)
    // are set to NULL, then only X and Y (or only X) are written.
    //   Shift corresponds to the begining of arrays' enumeration
    // X[Shift..Shift+len-1].
    int  i,l;
      l = Shift+len;
      for (i=Shift;i<l;i++)  {
        Write ( pstr("   ")  );
        Write ( X[i],MLength );
        if (Y)  {
          Write ( pstr(",   ") );
          Write ( Y[i],MLength );
        }
        if (Z)  {
          Write ( pstr(",   ") );
          Write ( Z[i],MLength );
        }
        LF();
      }
      return LF();
    }

    bool  File::WriteColumns ( rvector X, rvector Y,
                                   int len, int Shift, int MLength )  {
      return WriteColumns ( X,Y,NULL,len,Shift,MLength );
    }

    int  File::ReadColumns ( int maxlen, rvector X, rvector Y, rvector Z,
                              int xCol, int yCol, int zCol, int Shift )  {
    //   ReadColumns reads data stored by WriteColumns. X, Y, and Z must
    // be allocated prior to call.
    //   xCol, yCol and zCol specify the order number of columns
    // (starting from 0) to be read into X, Y and Z, correspondingly.
    // If zCol (or zCol and yCol) < 0 then Z (or Z and Y) are not read.
    //   Shift corresponds to the begining of arrays' enumeration
    // X[Shift..Shift+len-1].
    //   Returns number of lines read.
    int   DataLen;
    char  S[1025];
      DataLen = maxlen;
      _ReadColumns ( DataLen,S,sizeof(S),X,Y,Z,xCol,yCol,zCol,Shift );
      return DataLen;
    }

    int  File::ReadColumns ( int maxlen, rvector X, rvector Y,
                              int xCol, int yCol, int Shift )  {
      return  ReadColumns ( maxlen,X,Y,NULL,xCol,yCol,-1,Shift );
    }


    int  File::CreateReadColumns ( rvector & X, rvector & Y, rvector & Z,
                                    int xCol, int yCol, int zCol,
                                    int Shift )  {
    //   ReadColumns reads data stored by WriteColumns. X, Y, and Z
    // must be set to NULL prior to call. They will be allocated
    // within the procedure.
    //   xCol, yCol and zCol specify the order number of columns
    // (starting from 0) to be read into X, Y and Z, correspondingly.
    // If zCol (or zCol and yCol) < 0 then Z (or Z and Y) are not read.
    //   Shift corresponds to the begining of arrays' enumeration
    // X[Shift..Shift+len-1].
    //   Returns number of lines read, errors are reported by ErrorCode().
    int      i,j,DataLen;
    char     S[1025];
    bool  Ok;
      DataLen = 0;
      ErrCode = 0;
      if (!FileEnd())  {
        i = 0;
        j = 1;
        // the loop terminates at first blank line
        while ((i<j) && (!FileEnd()))  {
          j = ReadLine ( S,sizeof(S) );
          i = 0;
          // check for blank line
          while ((i<j) && (S[i]==' '))  i++;
          DataLen++;
        }
        if (i>=j)  DataLen--;
        if (DataLen>0)  {
          Ok = GetVectorMemory(X,DataLen,Shift);
          if (Ok && (yCol>=0))
            Ok = Ok && GetVectorMemory(Y,DataLen,Shift);
          if (Ok && (zCol>=0))
            Ok = Ok && GetVectorMemory(Z,DataLen,Shift);
          if (Ok)  {
            reset();
            _ReadColumns ( DataLen,S,sizeof(S),X,Y,Z,xCol,yCol,
                           zCol,Shift );
          } else  ErrCode = FileError_NoMemory;
        } else  ErrCode = FileError_NoDataFound;
      } else  ErrCode = FileError_NoDataFound;
      return DataLen;
    }

    int  File::CreateReadColumns ( rvector & X, rvector & Y,
                                    int xCol, int yCol, int Shift )  {
      return  CreateReadColumns ( X,Y,X,xCol,yCol,-1,Shift );
    }

    void  File::_ReadColumns ( int & DLen, pstr S, int SLen,
                                rvector X, rvector Y, rvector Z,
                                int xCol, int yCol, int zCol,
                                int Shift )  {
    int      i,is,j,k,m,n,cmax;
    char     SV[256];
    realtype Res;
      ErrCode = 0;
      i       = 0;
      cmax    = IMax(zCol,IMax(xCol,yCol));
      while ((i<DLen) && (ErrCode==0))  {
        k = ReadLine ( S,SLen );
        RemoveDelimiters(S,k);
        j = 0;
        m = -1;
        n = 0;
        while ((m<cmax) && (ErrCode==0))  {
          do {
            PickOutNumber ( S,SV,k,j );
            if ((m<0) && (SV[0]==char(0)) && (j>=k))  {
              DLen = i;
              return;
            }
            m++;
          } while ((m!=xCol) && (m!=yCol) && (m!=zCol));
          if (SV[0]==char(0))  {
            if (n>0)  ErrCode = FileError_NoColumn;
                else  ErrCode = FileError_ShortData;
          } else  {
            Res = GetNumber ( SV );
            if (ErrCode==0)  {
              is = i+Shift;
              if      (m==xCol)  X[is] = Res;
              else if (m==yCol)  Y[is] = Res;
                           else  Z[is] = Res;
              n++;
            }
          }
        }
        if ((ErrCode==0) && (n<2))  ErrCode = FileError_NoColumn;
        i++;
      }
      if ((ErrCode==FileError_ShortData) && (i>1))  {
        ErrCode = 0;
        DLen    = i-1;
      }
      if (ErrCode!=0)  ErrCode = FileError_BadData;
    }

    void  RemoveDelimiters ( pstr S, int SLen )  {
    int j;
      for (j=0;j<SLen;j++)
        if ((S[j]==',') || (S[j]==';') ||
            (S[j]==':') || (S[j]==char(9)))
          S[j] = ' ';
    }

    void  PickOutNumber ( cpstr S, pstr SV, int SLen, int & j )  {
    int l;
      l = 0;
      while ((j<SLen) && (S[j]==' '))  j++;
      if ((S[j]=='+') || (S[j]=='-'))  SV[l++] = S[j++];
      if (S[j]=='.')                   SV[l++] = '0';
      while ((j<SLen) && (S[j]!=' '))  SV[l++] = S[j++];
      SV[l] = char(0);
    }

    realtype  File::GetNumber ( cpstr S )  {
    char *   endptr;
    realtype V;
      V = strtod ( S,&endptr );
      if ((*endptr!=' ') && (*endptr))
        ErrCode = FileError_BadData;
      return V;
    }

    cpstr FileError ( int ErrCode )  {
      switch (ErrCode)  {
        case 0                      : return "Ok";
        case FileError_NoMemory     : return "Insufficient memory";
        case FileError_NoDataFound  : return "No data found";
        case FileError_NoColumn     : return "No column structure";
        case FileError_BadData      : return "Incorrect data format";
        case FileError_WrongMemoryAllocation
                                    : return "Wrong Memory Allocation";
        default : return "Unknown I/O error";
      }
    }

 }

}
