/* spacegroup_data.cpp: spacegroup data tables */
//C Copyright (C) 2000-2006 Kevin Cowtan and University of York
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA
/* Some code coverted from the paper of Hall & Grosse-Kunstleve */


#include "spacegroup_data.h"


namespace clipper {


namespace data {


  // ASU fns
  bool ASU_111( const int& h, const int& k, const int& l )
    { return (l>0 || (l==0 && (h>0 || (h==0 && k>=0)))); }
  bool ASU_112( const int& h, const int& k, const int& l )
    { return (l>=0 && (h>0 || (h==0 && k>=0))); }
  bool ASU_121( const int& h, const int& k, const int& l )
    { return (k>=0 && (l>0 || (l==0 && h>=0))); }
  bool ASU_211( const int& h, const int& k, const int& l )
    { return (h>=0 && (k>0 || (k==0 && l>=0))); }
  bool ASU_21U( const int& h, const int& k, const int& l )
    { return (h+k>=0 && (l>0 || (l==0 && h-k>=0))); }
  bool ASU_21V( const int& h, const int& k, const int& l )
    { return (l+h>=0 && (k>0 || (k==0 && l-h>=0))); }
  bool ASU_21W( const int& h, const int& k, const int& l )
    { return (k+l>=0 && (h>0 || (h==0 && k-l>=0))); }
  bool ASU_21X( const int& h, const int& k, const int& l )
    { return (h-k>=0 && (l>0 || (l==0 && h+k>=0))); }
  bool ASU_21Y( const int& h, const int& k, const int& l )
    { return (l-h>=0 && (k>0 || (k==0 && l+h>=0))); }
  bool ASU_21Z( const int& h, const int& k, const int& l )
    { return (k-l>=0 && (h>0 || (h==0 && k+l>=0))); }
  bool ASU_222( const int& h, const int& k, const int& l )
    { return (h>=0 && k>=0 && l>=0); }
  bool ASU_22U( const int& h, const int& k, const int& l )
    { return (h<=k && h>=-k && l>=0); }
  bool ASU_22V( const int& h, const int& k, const int& l )
    { return (l<=h && l>=-h && k>=0); }
  bool ASU_22W( const int& h, const int& k, const int& l )
    { return (k<=l && k>=-l && h>=0); }
  bool ASU_114( const int& h, const int& k, const int& l )
    { return (l>=0 && ((h>=0 && k>0) || (h==0 && k==0))); }
  bool ASU_141( const int& h, const int& k, const int& l )
    { return (k>=0 && ((l>=0 && h>0) || (l==0 && h==0))); }
  bool ASU_411( const int& h, const int& k, const int& l )
    { return (h>=0 && ((k>=0 && l>0) || (k==0 && l==0))); }
  bool ASU_224( const int& h, const int& k, const int& l )
    { return (h>=k && k>=0 && l>=0); }
  bool ASU_242( const int& h, const int& k, const int& l )
    { return (l>=h && h>=0 && k>=0); }
  bool ASU_422( const int& h, const int& k, const int& l )
    { return (k>=l && l>=0 && h>=0); }
  bool ASU_113( const int& h, const int& k, const int& l )
    { return (h>=0 && k>0) || (h==0 && k==0 && l >= 0); }
  bool ASU_131( const int& h, const int& k, const int& l )
    { return (l>=0 && h>0) || (l==0 && h==0 && k >= 0); }
  bool ASU_311( const int& h, const int& k, const int& l )
    { return (k>=0 && l>0) || (k==0 && l==0 && h >= 0); }
  bool ASU_11T( const int& h, const int& k, const int& l )
    { return (h<=0 && k>0) || (h==0 && k==0 && l >= 0); }
  bool ASU_1T1( const int& h, const int& k, const int& l )
    { return (l<=0 && h>0) || (l==0 && h==0 && k >= 0); }
  bool ASU_T11( const int& h, const int& k, const int& l )
    { return (k<=0 && l>0) || (k==0 && l==0 && h >= 0); }
  bool ASU_31A( const int& h, const int& k, const int& l )
    { return (k-l>=0 && l-h>0) || (h==l && k==l && h+k+l>=0); }
  bool ASU_31B( const int& h, const int& k, const int& l )
    { return (k-l>=0 && l+h>0) || (-h==l && k==l && -h+k+l>=0); }
  bool ASU_31C( const int& h, const int& k, const int& l )
    { return (-k-l>=0 && l-h>0) || (h==l && -k==l && h-k+l>=0); }
  bool ASU_31D( const int& h, const int& k, const int& l )
    { return (k+l>=0 && -l-h>0) || (h==-l && k==-l && h+k-l>=0); }
  bool ASU_223( const int& h, const int& k, const int& l )
    { return (h>=k && k>=0 && (k>0 || l>=0)); }
  bool ASU_232( const int& h, const int& k, const int& l )
    { return (l>=h && h>=0 && (h>0 || k>=0)); }
  bool ASU_322( const int& h, const int& k, const int& l )
    { return (k>=l && l>=0 && (l>0 || h>=0)); }
  bool ASU_32A( const int& h, const int& k, const int& l )
    { return (h>=k && k+l>=h+h && (k+l>h+h || h+k+l>=0)); }
  bool ASU_32B( const int& h, const int& k, const int& l )
    { return (-h>=k && k+l>=-h-h && (k+l>-h-h || -h+k+l>=0)); }
  bool ASU_32C( const int& h, const int& k, const int& l )
    { return (h>=-k && -k+l>=h+h && (-k+l>h+h || h-k+l>=0)); }
  bool ASU_32D( const int& h, const int& k, const int& l )
    { return (h>=k && k-l>=h+h && (k-l>h+h || h+k-l>=0)); }
  bool ASU_32U( const int& h, const int& k, const int& l )
    { return (h>=k && k>=0 && (h>k || l>=0)); }
  bool ASU_32V( const int& h, const int& k, const int& l )
    { return (k>=l && l>=0 && (k>l || h>=0)); }
  bool ASU_32W( const int& h, const int& k, const int& l )
    { return (l>=h && h>=0 && (l>h || k>=0)); }
  bool ASU_32X( const int& h, const int& k, const int& l )
    { return (-h>=k && k>=0 && (-h>k || l>=0)); }
  bool ASU_32Y( const int& h, const int& k, const int& l )
    { return (-k>=l && l>=0 && (-k>l || h>=0)); }
  bool ASU_32Z( const int& h, const int& k, const int& l )
    { return (-l>=h && h>=0 && (-l>h || k>=0)); }
  bool ASU_M3B( const int& h, const int& k, const int& l )
    { return (h>=0 && ((l>=h && k>h) || (l==h && k==h))); }	   
  bool ASU_M3M( const int& h, const int& k, const int& l )
    { return (k>=l && l>=h && h>=0); }

// ABCD are 4 dirns for body diagonal, UVWXYZ are 6 dirns for face diagonal
LGdata lgdata[] = {
  { 0x2df60a45, &ASU_111, "-1" },     // -P 1		   
  { 0x9b92779f, &ASU_112, "2/m" },    // -P 2		   
  { 0x068dfeb3, &ASU_121, "2/m" },    // -P 2y		   
  { 0xb08d46cf, &ASU_211, "2/m" },    // -P 2x		   
  { 0x12339040, &ASU_21U, "2/m" },    // -P 2"		   
  { 0x44aa9a14, &ASU_21V, "2/m" },    // -P 2y"	   
  { 0x53e4b366, &ASU_21W, "2/m" },    // -P 2x"	   
  { 0x4321a07d, &ASU_21X, "2/m" },    // -P 2'		   
  { 0x1c1b5411, &ASU_21Y, "2/m" },    // -P 2y'	   
  { 0xe34a99ed, &ASU_21Z, "2/m" },    // -P 2x'	   
  { 0xe7243bbc, &ASU_222, "mmm" },    // -P 2 2	   
  { 0x885920bf, &ASU_22U, "mmm" },    // -P 2 2"	   
  { 0xe980f874, &ASU_22V, "mmm" },    // -P 2 2" (y,z,x)	   
  { 0x05c7f86e, &ASU_22W, "mmm" },    // -P 2 2" (z,x,y)	   
  { 0xfdd759b5, &ASU_113, "-3" },     // -P 3		   
  { 0xf2769c28, &ASU_131, "-3" },     // -P 3 (y,z,x)	   
  { 0xcd4b8428, &ASU_311, "-3" },     // -P 3 (z,x,y)	   
  { 0x07fa5ca1, &ASU_11T, "-3" },     // -P 3 (-x,y,z)
  { 0x0f070468, &ASU_1T1, "-3" },     // -P 3 (y,z,-x)
  { 0x9bd3dcf4, &ASU_T11, "-3" },     // -P 3 (z,-x,y)
  { 0xd9a29bac, &ASU_31A, "-3" },     // -P 3*		   
  { 0x627a2f8c, &ASU_31B, "-3" },     // -P 3* (-x,y,z)   
  { 0xcfb71d71, &ASU_31C, "-3" },     // -P 3* (x,-y,z)   
  { 0x8a86426e, &ASU_31D, "-3" },     // -P 3* (x,y,-z)   
  { 0xf74c7f83, &ASU_223, "-3m" },    // -P 3 2	   
  { 0x573b981c, &ASU_232, "-3m" },    // -P 3 2 (y,z,x)   
  { 0x1799544d, &ASU_322, "-3m" },    // -P 3 2 (z,x,y)   
  { 0x1c80e47a, &ASU_32A, "-3m" },    // -P 3* 2	   
  { 0xea7284da, &ASU_32B, "-3m" },    // -P 3* 2 (-x,y,z) 
  { 0xb193db73, &ASU_32C, "-3m" },    // -P 3* 2 (x,-y,z) 
  { 0x04fecdf9, &ASU_32D, "-3m" },    // -P 3* 2 (-x,-y,z)
  { 0xfc3edafb, &ASU_32U, "-3m" },    // -P 3 2"	   
  { 0xd60d11a0, &ASU_32V, "-3m" },    // -P 3 2" (z,x,y)  
  { 0xf7d5112f, &ASU_32W, "-3m" },    // -P 3 2" (y,z,x)  
  { 0xfbb8f18b, &ASU_32X, "-3m" },    // -P 3 2" (-x,y,z) 
  { 0x530fcba9, &ASU_32Y, "-3m" },    // -P 3 2" (z,-x,y) 
  { 0xa3d49592, &ASU_32Z, "-3m" },    // -P 3 2" (y,z,-x) 
  { 0x7f473453, &ASU_114, "4/m" },    // -P 4		   
  { 0x081d78e5, &ASU_141, "4/m" },    // -P 4 (y,z,x)	   
  { 0x6a9fa6e5, &ASU_411, "4/m" },    // -P 4 (z,x,y)	   
  { 0xb8f2113e, &ASU_224, "4/mmm" },  // -P 4 2	   
  { 0xb9ce4369, &ASU_242, "4/mmm" },  // -P 4 2 (y,z,x)   
  { 0xc39787b7, &ASU_422, "4/mmm" },  // -P 4 2 (z,x,y)   
  { 0x32dacfb6, &ASU_114, "6/m" },    // -P 6		   
  { 0xe7987b0c, &ASU_141, "6/m" },    // -P 6 (y,z,x)	   
  { 0xb5b69658, &ASU_411, "6/m" },    // -P 6 (z,x,y)	   
  { 0xf1fc7952, &ASU_224, "6/mmm" },  // -P 6 2	   
  { 0x386f0ab4, &ASU_242, "6/mmm" },  // -P 6 2 (y,z,x)   
  { 0xcd531b66, &ASU_422, "6/mmm" },  // -P 6 2 (z,x,y)   
  { 0x72e55913, &ASU_M3B, "m-3" },    // -P 2 2 3	   
  { 0x74c407d3, &ASU_M3M, "m-3m" }    // -P 4 2 3         
};
int lgdata_size = sizeof( lgdata ) / sizeof( lgdata[0] );

// From Hall & Grosse-Kunstleve, 'Concise Space group symbols' 
SGdata sgdata[] = {
  {0xc704dd7b,"P 1"             ,"P 1"       ,' ',  1},
  {0x2df60a45,"-P 1"            ,"P -1"      ,' ',  2},
  {0x902cf668,"P 2y"            ,"P 1 2 1"   ,' ',  3},
  {0x90a34743,"P 2"             ,"P 1 1 2"   ,' ',  3},
  {0x9b02bfbf,"P 2x"            ,"P 2 1 1"   ,' ',  3},
  {0xda168154,"P 2yb"           ,"P 1 21 1"  ,' ',  4},
  {0x7851c0d1,"P 2c"            ,"P 1 1 21"  ,' ',  4},
  {0xae0e24db,"P 2xa"           ,"P 21 1 1"  ,' ',  4},
  {0xdeb5173b,"C 2y"            ,"C 1 2 1"   ,' ',  5},
  {0x1e4dbca6,"A 2y"            ,"A 1 2 1"   ,' ',  5},
  {0x6042abb6,"I 2y"            ,"I 1 2 1"   ,' ',  5},
  {0x2a55a450,"A 2"             ,"A 1 1 2"   ,' ',  5},
  {0x458c0b8c,"B 2"             ,"B 1 1 2"   ,' ',  5},
  {0x545ab340,"I 2"             ,"I 1 1 2"   ,' ',  5},
  {0xc4ef9a2d,"B 2x"            ,"B 2 1 1"   ,' ',  5},
  {0x6bce9e6c,"C 2x"            ,"C 2 1 1"   ,' ',  5},
  {0xd53922e1,"I 2x"            ,"I 2 1 1"   ,' ',  5},
  {0x93121e4d,"P -2y"           ,"P 1 m 1"   ,' ',  6},
  {0x939daf66,"P -2"            ,"P 1 1 m"   ,' ',  6},
  {0x6dd7504b,"P -2x"           ,"P m 1 1"   ,' ',  6},
  {0x7be099df,"P -2yc"          ,"P 1 c 1"   ,' ',  7},
  {0x4eec02bb,"P -2yac"         ,"P 1 n 1"   ,' ',  7},
  {0xa61e8529,"P -2ya"          ,"P 1 a 1"   ,' ',  7},
  {0xa6913402,"P -2a"           ,"P 1 1 a"   ,' ',  7},
  {0xecab433e,"P -2ab"          ,"P 1 1 n"   ,' ',  7},
  {0xd9a7d85a,"P -2b"           ,"P 1 1 b"   ,' ',  7},
  {0x27ed2777,"P -2xb"          ,"P b 1 1"   ,' ',  7},
  {0x2ac91f43,"P -2xbc"         ,"P n 1 1"   ,' ',  7},
  {0x8525d7d9,"P -2xc"          ,"P c 1 1"   ,' ',  7},
  {0x77b286e1,"C -2y"           ,"C 1 m 1"   ,' ',  8},
  {0xb74a2d7c,"A -2y"           ,"A 1 m 1"   ,' ',  8},
  {0xc9453a6c,"I -2y"           ,"I 1 m 1"   ,' ',  8},
  {0x8352358a,"A -2"            ,"A 1 1 m"   ,' ',  8},
  {0xec8b9a56,"B -2"            ,"B 1 1 m"   ,' ',  8},
  {0xfd5d229a,"I -2"            ,"I 1 1 m"   ,' ',  8},
  {0xddc257bb,"B -2x"           ,"B m 1 1"   ,' ',  8},
  {0x72e353fa,"C -2x"           ,"C m 1 1"   ,' ',  8},
  {0xcc14ef77,"I -2x"           ,"I m 1 1"   ,' ',  8},
  {0x38ef08aa,"C -2yc"          ,"C 1 c 1"   ,' ',  9},
  {0x54faf6d0,"A -2yac"         ,"A 1 n 1"   ,' ',  9},
  {0x9d54258d,"I -2ya"          ,"I 1 a 1"   ,' ',  9},
  {0xe35b329d,"A -2ya"          ,"A 1 a 1"   ,' ',  9},
  {0x6cfe174b,"C -2ybc"         ,"C 1 n 1"   ,' ',  9},
  {0x2af5e1c0,"I -2yc"          ,"I 1 c 1"   ,' ',  9},
  {0xd7432a6b,"A -2a"           ,"A 1 1 a"   ,' ',  9},
  {0x0f3b41fa,"B -2bc"          ,"B 1 1 n"   ,' ',  9},
  {0x4afce6d7,"I -2b"           ,"I 1 1 b"   ,' ',  9},
  {0x5b2a5e1b,"B -2b"           ,"B 1 1 b"   ,' ',  9},
  {0x60e2ee26,"A -2ac"          ,"A 1 1 n"   ,' ',  9},
  {0xa94c3d7b,"I -2a"           ,"I 1 1 a"   ,' ',  9},
  {0x6a6393f6,"B -2xb"          ,"B b 1 1"   ,' ',  9},
  {0x69afc250,"C -2xbc"         ,"C n 1 1"   ,' ',  9},
  {0x2fa434db,"I -2xc"          ,"I c 1 1"   ,' ',  9},
  {0x3dbeddb1,"C -2xc"          ,"C c 1 1"   ,' ',  9},
  {0x3e728c17,"B -2xbc"         ,"B n 1 1"   ,' ',  9},
  {0x7bb52b3a,"I -2xb"          ,"I b 1 1"   ,' ',  9},
  {0x068dfeb3,"-P 2y"           ,"P 1 2/m 1" ,' ', 10},
  {0x9b92779f,"-P 2"            ,"P 1 1 2/m" ,' ', 10},
  {0xb08d46cf,"-P 2x"           ,"P 2/m 1 1" ,' ', 10},
  {0x54fa8558,"-P 2yb"          ,"P 1 21/m 1",' ', 11},
  {0x31194672,"-P 2c"           ,"P 1 1 21/m",' ', 11},
  {0xce8251df,"-P 2xa"          ,"P 21/m 1 1",' ', 11},
  {0x09efdd05,"-C 2y"           ,"C 1 2/m 1" ,' ', 12},
  {0x8e118804,"-A 2y"           ,"A 1 2/m 1" ,' ', 12},
  {0x5c5d4c9f,"-I 2y"           ,"I 1 2/m 1" ,' ', 12},
  {0x5a48cc10,"-A 2"            ,"A 1 1 2/m" ,' ', 12},
  {0x97e84b5c,"-B 2"            ,"B 1 1 2/m" ,' ', 12},
  {0x8804088b,"-I 2"            ,"I 1 1 2/m" ,' ', 12},
  {0x05a3ecad,"-B 2x"           ,"B 2/m 1 1" ,' ', 12},
  {0x4ffd3ee0,"-C 2x"           ,"C 2/m 1 1" ,' ', 12},
  {0x1a4faf7a,"-I 2x"           ,"I 2/m 1 1" ,' ', 12},
  {0xac06cf5e,"-P 2yc"          ,"P 1 2/c 1" ,' ', 13},
  {0xf817d0bf,"-P 2yac"         ,"P 1 2/n 1" ,' ', 13},
  {0x529ce152,"-P 2ya"          ,"P 1 2/a 1" ,' ', 13},
  {0xcf83687e,"-P 2a"           ,"P 1 1 2/a" ,' ', 13},
  {0x9df41395,"-P 2ab"          ,"P 1 1 2/n" ,' ', 13},
  {0xc9e50c74,"-P 2b"           ,"P 1 1 2/b" ,' ', 13},
  {0x4aada62e,"-P 2xb"          ,"P 2/b 1 1" ,' ', 13},
  {0xf45a1aa3,"-P 2xbc"         ,"P 2/n 1 1" ,' ', 13},
  {0xe58ca26f,"-P 2xc"          ,"P 2/c 1 1" ,' ', 13},
  {0x41572403,"-P 2ybc"         ,"P 1 21/c 1",' ', 14},
  {0x15463be2,"-P 2yn"          ,"P 1 21/n 1",' ', 14},
  {0x00eb9ab9,"-P 2yab"         ,"P 1 21/a 1",' ', 14},
  {0x65085993,"-P 2ac"          ,"P 1 1 21/a",' ', 14},
  {0x8859b2ce,"-P 2n"           ,"P 1 1 21/n",' ', 14},
  {0xdc48ad2f,"-P 2bc"          ,"P 1 1 21/b",' ', 14},
  {0x34a2b13e,"-P 2xab"         ,"P 21/b 1 1",' ', 14},
  {0x8a550db3,"-P 2xn"          ,"P 21/n 1 1",' ', 14},
  {0x9b83b57f,"-P 2xac"         ,"P 21/c 1 1",' ', 14},
  {0x673d531e,"-C 2yc"          ,"C 1 2/c 1" ,' ', 15},
  {0x428e62a8,"-A 2yac"         ,"A 1 2/n 1" ,' ', 15},
  {0xb7f841ea,"-I 2ya"          ,"I 1 2/a 1" ,' ', 15},
  {0x65b48571,"-A 2ya"          ,"A 1 2/a 1" ,' ', 15},
  {0x8c985e6b,"-C 2ybc"         ,"C 1 2/n 1" ,' ', 15},
  {0x90c2a633,"-I 2yc"          ,"I 1 2/c 1" ,' ', 15},
  {0xb1edc165,"-A 2a"           ,"A 1 1 2/a" ,' ', 15},
  {0x5b77a1f0,"-B 2bc"          ,"B 1 1 2/n" ,' ', 15},
  {0xaf3eef52,"-I 2b"           ,"I 1 1 2/b" ,' ', 15},
  {0xb0d2ac85,"-B 2b"           ,"B 1 1 2/b" ,' ', 15},
  {0x96d726bc,"-A 2ac"          ,"A 1 1 2/n" ,' ', 15},
  {0x63a105fe,"-I 2a"           ,"I 1 1 2/a" ,' ', 15},
  {0xe4301506,"-B 2xb"          ,"B 2/b 1 1" ,' ', 15},
  {0xaecdc717,"-C 2xbc"         ,"C 2/n 1 1" ,' ', 15},
  {0x62f43297,"-I 2xc"          ,"I 2/c 1 1" ,' ', 15},
  {0x37e5a351,"-C 2xc"          ,"C 2/c 1 1" ,' ', 15},
  {0x7d187140,"-B 2xbc"         ,"B 2/n 1 1" ,' ', 15},
  {0xfbdc56d1,"-I 2xb"          ,"I 2/b 1 1" ,' ', 15},
  {0xa959fc0b,"P 2 2"           ,"P 2 2 2"   ,' ', 16},
  {0x03d2cde6,"P 2c 2"          ,"P 2 2 21"  ,' ', 17},
  {0xd756eb1b,"P 2a 2a"         ,"P 21 2 2"  ,' ', 17},
  {0x010e6701,"P 2 2b"          ,"P 2 21 2"  ,' ', 17},
  {0x2b106ff0,"P 2 2ab"         ,"P 21 21 2" ,' ', 18},
  {0xee8326bb,"P 2bc 2"         ,"P 2 21 21" ,' ', 18},
  {0x82570fbb,"P 2ac 2ac"       ,"P 21 2 21" ,' ', 18},
  {0x8f7a6eec,"P 2ac 2ab"       ,"P 21 21 21",' ', 19},
  {0xd3dceae0,"C 2c 2"          ,"C 2 2 21"  ,' ', 20},
  {0xa3d855bc,"A 2a 2a"         ,"A 21 2 2"  ,' ', 20},
  {0x31f9a8c4,"B 2 2b"          ,"B 2 21 2"  ,' ', 20},
  {0xbd0e64fb,"C 2 2"           ,"C 2 2 2"   ,' ', 21},
  {0x3af031fa,"A 2 2"           ,"A 2 2 2"   ,' ', 21},
  {0xf750b6b6,"B 2 2"           ,"B 2 2 2"   ,' ', 21},
  {0x8a25cd68,"F 2 2"           ,"F 2 2 2"   ,' ', 22},
  {0xe8bcf561,"I 2 2"           ,"I 2 2 2"   ,' ', 23},
  {0x7ba265f9,"I 2b 2c"         ,"I 21 21 21",' ', 24},
  {0x715a8a38,"P 2 -2"          ,"P m m 2"   ,' ', 25},
  {0x34467527,"P -2 2"          ,"P 2 m m"   ,' ', 25},
  {0x454292ce,"P -2 -2"         ,"P m 2 m"   ,' ', 25},
  {0xdbd1bbd5,"P 2c -2"         ,"P m c 21"  ,' ', 26},
  {0x8ed05f75,"P 2c -2c"        ,"P c m 21"  ,' ', 26},
  {0x1e587dd6,"P -2a 2a"        ,"P 21 m a"  ,' ', 26},
  {0x4a496237,"P -2 2a"         ,"P 21 a m"  ,' ', 26},
  {0xed1509c4,"P -2 -2b"        ,"P b 21 m"  ,' ', 26},
  {0x1735e925,"P -2b -2"        ,"P m 21 b"  ,' ', 26},
  {0x245b6e98,"P 2 -2c"         ,"P c c 2"   ,' ', 27},
  {0x60576ac6,"P -2a 2"         ,"P 2 a a"   ,' ', 27},
  {0xbf62722f,"P -2b -2b"       ,"P b 2 b"   ,' ', 27},
  {0x0f559d28,"P 2 -2a"         ,"P m a 2"   ,' ', 28},
  {0x8b7a6ad9,"P 2 -2b"         ,"P b m 2"   ,' ', 28},
  {0x66310ecc,"P -2b 2"         ,"P 2 m b"   ,' ', 28},
  {0x9ecd44ca,"P -2c 2"         ,"P 2 c m"   ,' ', 28},
  {0x1043766e,"P -2c -2c"       ,"P c 2 m"   ,' ', 28},
  {0x3b4d85de,"P -2a -2a"       ,"P m 2 a"   ,' ', 28},
  {0xf0df4865,"P 2c -2ac"       ,"P c a 21"  ,' ', 29},
  {0xc427e492,"P 2c -2b"        ,"P b c 21"  ,' ', 29},
  {0x183e19dc,"P -2b 2a"        ,"P 21 a b"  ,' ', 29},
  {0xb4d34c3b,"P -2ac 2a"       ,"P 21 c a"  ,' ', 29},
  {0xa7e2b223,"P -2bc -2c"      ,"P c 21 b"  ,' ', 29},
  {0x931a1ed4,"P -2a -2ab"      ,"P b 21 a"  ,' ', 29},
  {0x358dd654,"P 2 -2bc"        ,"P n c 2"   ,' ', 30},
  {0x5a547988,"P 2 -2ac"        ,"P c n 2"   ,' ', 30},
  {0xcadc5b2b,"P -2ac 2"        ,"P 2 n a"   ,' ', 30},
  {0x3220112d,"P -2ab 2"        ,"P 2 a n"   ,' ', 30},
  {0xc16d653f,"P -2ab -2ab"     ,"P b 2 n"   ,' ', 30},
  {0x0195cea2,"P -2bc -2bc"     ,"P n 2 b"   ,' ', 30},
  {0x8fc0a434,"P 2ac -2"        ,"P m n 21"  ,' ', 31},
  {0x72570ce4,"P 2bc -2bc"      ,"P n m 21"  ,' ', 31},
  {0xb60fe6dc,"P -2ab 2ab"      ,"P 21 m n"  ,' ', 31},
  {0x1f488697,"P -2 2ac"        ,"P 21 n m"  ,' ', 31},
  {0x464f1412,"P -2 -2bc"       ,"P n 21 m"  ,' ', 31},
  {0x4324f6c4,"P -2ab -2"       ,"P m 21 n"  ,' ', 31},
  {0xf5757dc9,"P 2 -2ab"        ,"P b a 2"   ,' ', 32},
  {0x739caf97,"P -2bc 2"        ,"P 2 c b"   ,' ', 32},
  {0x6e4c617e,"P -2ac -2ac"     ,"P c 2 a"   ,' ', 32},
  {0x04df4f0f,"P 2c -2n"        ,"P n a 21"  ,' ', 33},
  {0xba28f382,"P 2c -2ab"       ,"P b n 21"  ,' ', 33},
  {0x0d93b887,"P -2bc 2a"       ,"P 21 n b"  ,' ', 33},
  {0x5982a766,"P -2n 2a"        ,"P 21 c n"  ,' ', 33},
  {0xd9eda533,"P -2n -2ac"      ,"P c 21 n"  ,' ', 33},
  {0xc83b1dff,"P -2ac -2n"      ,"P n 21 a"  ,' ', 33},
  {0x4b82c144,"P 2 -2n"         ,"P n n 2"   ,' ', 34},
  {0x278db076,"P -2n 2"         ,"P 2 n n"   ,' ', 34},
  {0x7f9ad9b2,"P -2n -2n"       ,"P n 2 n"   ,' ', 34},
  {0x93d6f19f,"C 2 -2"          ,"C m m 2"   ,' ', 35},
  {0xeea975ee,"A -2 2"          ,"A 2 m m"   ,' ', 35},
  {0xb633cb54,"B -2 -2"         ,"B m 2 m"   ,' ', 35},
  {0xfd047f84,"C 2c -2"         ,"C m c 21"  ,' ', 36},
  {0x851ce235,"C 2c -2c"        ,"C c m 21"  ,' ', 36},
  {0x9c241cdd,"A -2a 2a"        ,"A 21 m a"  ,' ', 36},
  {0x778111a8,"A -2 2a"         ,"A 21 a m"  ,' ', 36},
  {0x709ad526,"B -2 -2b"        ,"B b 21 m"  ,' ', 36},
  {0x91092c8d,"B -2b -2"        ,"B m 21 b"  ,' ', 36},
  {0xebce6c2e,"C 2 -2c"         ,"C c c 2"   ,' ', 37},
  {0x050c789b,"A -2a 2"         ,"A 2 a a"   ,' ', 37},
  {0x57a032ff,"B -2b -2b"       ,"B b 2 b"   ,' ', 37},
  {0x1428a49e,"A 2 -2"          ,"A m m 2"   ,' ', 38},
  {0xd98823d2,"B 2 -2"          ,"B m m 2"   ,' ', 38},
  {0x2309f2a2,"B -2 2"          ,"B 2 m m"   ,' ', 38},
  {0x695720ef,"C -2 2"          ,"C 2 m m"   ,' ', 38},
  {0xfc6d1919,"C -2 -2"         ,"C m 2 m"   ,' ', 38},
  {0x7b934c18,"A -2 -2"         ,"A m 2 m"   ,' ', 38},
  {0xf5bb5d35,"A 2 -2c"         ,"A b m 2"   ,' ', 39},
  {0x40a04794,"B 2 -2c"         ,"B m a 2"   ,' ', 39},
  {0xc8acffd7,"B -2c 2"         ,"B 2 c m"   ,' ', 39},
  {0x82f22d9a,"C -2b 2"         ,"C 2 m b"   ,' ', 39},
  {0x65457d5f,"C -2b -2b"       ,"C m 2 a"   ,' ', 39},
  {0x9a00b5b3,"A -2c -2c"       ,"A c 2 m"   ,' ', 39},
  {0x8d00c0d8,"A 2 -2a"         ,"A m a 2"   ,' ', 40},
  {0x381bda79,"B 2 -2b"         ,"B b m 2"   ,' ', 40},
  {0x0433157b,"B -2b 2"         ,"B 2 m b"   ,' ', 40},
  {0x0785aef4,"C -2c 2"         ,"C 2 c m"   ,' ', 40},
  {0x847584a8,"C -2c -2c"       ,"C c 2 m"   ,' ', 40},
  {0xe2bb285e,"A -2a -2a"       ,"A m 2 a"   ,' ', 40},
  {0x6c933973,"A 2 -2ac"        ,"A b a 2"   ,' ', 41},
  {0xa133be3f,"B 2 -2bc"        ,"B b a 2"   ,' ', 41},
  {0xef96180e,"B -2bc 2"        ,"B 2 c b"   ,' ', 41},
  {0xec20a381,"C -2bc 2"        ,"C 2 c b"   ,' ', 41},
  {0x1d5de0ee,"C -2bc -2bc"     ,"C c 2 a"   ,' ', 41},
  {0x0328d1f5,"A -2ac -2ac"     ,"A c 2 a"   ,' ', 41},
  {0x18172b53,"F 2 -2"          ,"F m m 2"   ,' ', 42},
  {0xad1a17db,"F -2 2"          ,"F 2 m m"   ,' ', 42},
  {0xa88fc3e9,"F -2 -2"         ,"F m 2 m"   ,' ', 42},
  {0x921865a5,"F 2 -2d"         ,"F d d 2"   ,' ', 43},
  {0x1b59dbd3,"F -2d 2"         ,"F 2 d d"   ,' ', 43},
  {0x22808d1f,"F -2d -2d"       ,"F d 2 d"   ,' ', 43},
  {0xc6646005,"I 2 -2"          ,"I m m 2"   ,' ', 44},
  {0x3ce5b175,"I -2 2"          ,"I 2 m m"   ,' ', 44},
  {0xa9df8883,"I -2 -2"         ,"I m 2 m"   ,' ', 44},
  {0xbedffde8,"I 2 -2c"         ,"I b a 2"   ,' ', 45},
  {0xd740bc00,"I -2a 2"         ,"I 2 c b"   ,' ', 45},
  {0x484c7128,"I -2b -2b"       ,"I c 2 a"   ,' ', 45},
  {0x5f4c0443,"I 2 -2a"         ,"I m a 2"   ,' ', 46},
  {0x27f799ae,"I 2 -2b"         ,"I b m 2"   ,' ', 46},
  {0x1bdf56ac,"I -2b 2"         ,"I 2 m b"   ,' ', 46},
  {0xf07a5bd9,"I -2c 2"         ,"I 2 c m"   ,' ', 46},
  {0xd164156e,"I -2c -2c"       ,"I c 2 m"   ,' ', 46},
  {0x30f7ecc5,"I -2a -2a"       ,"I m 2 a"   ,' ', 46},
  {0xe7243bbc,"-P 2 2"          ,"P m m m"   ,' ', 47},
  {0x51e85b6c,"P 2 2 -1n"       ,"P n n n"   ,'1', 48},
  {0x78ff5b81,"-P 2ab 2bc"      ,"P n n n"   ,'2', 48},
  {0x353736e5,"-P 2 2c"         ,"P c c m"   ,' ', 49},
  {0x0c8136c9,"-P 2a 2"         ,"P m a a"   ,' ', 49},
  {0x02b17898,"-P 2b 2b"        ,"P b m b"   ,' ', 49},
  {0x1392524e,"P 2 2 -1ab"      ,"P b a n"   ,'1', 50},
  {0xe91475ed,"-P 2ab 2b"       ,"P b a n"   ,'2', 50},
  {0x4e9de450,"P 2 2 -1bc"      ,"P n c b"   ,'1', 50},
  {0x935a56f4,"-P 2b 2bc"       ,"P n c b"   ,'2', 50},
  {0x38bbae1f,"P 2 2 -1ac"      ,"P c n a"   ,'1', 50},
  {0xde923b90,"-P 2a 2c"        ,"P c n a"   ,'2', 50},
  {0x6179e0c6,"-P 2a 2a"        ,"P m m a"   ,' ', 51},
  {0x2be88448,"-P 2b 2"         ,"P m m b"   ,' ', 51},
  {0xce7dc76c,"-P 2 2b"         ,"P b m m"   ,' ', 51},
  {0xb013e0d3,"-P 2c 2c"        ,"P c m m"   ,' ', 51},
  {0x6200ed8a,"-P 2c 2"         ,"P m c m"   ,' ', 51},
  {0x8adcedb3,"-P 2 2a"         ,"P m a m"   ,' ', 51},
  {0xc514e0e2,"-P 2a 2bc"       ,"P n n a"   ,' ', 52},
  {0xfea280fb,"-P 2b 2n"        ,"P n n b"   ,' ', 52},
  {0x38e06b40,"-P 2n 2b"        ,"P b n n"   ,' ', 52},
  {0x637980f3,"-P 2ab 2c"       ,"P c n n"   ,' ', 52},
  {0x15078d8e,"-P 2ab 2n"       ,"P n c n"   ,' ', 52},
  {0xa90b452c,"-P 2n 2bc"       ,"P n a n"   ,' ', 52},
  {0x89a5e0ff,"-P 2ac 2"        ,"P m n a"   ,' ', 53},
  {0x42ae4859,"-P 2bc 2bc"      ,"P n m b"   ,' ', 53},
  {0x84eca3e2,"-P 2ab 2ab"      ,"P b m n"   ,' ', 53},
  {0x58cfe0ea,"-P 2 2ac"        ,"P c n m"   ,' ', 53},
  {0x2eb1ed97,"-P 2 2bc"        ,"P n c m"   ,' ', 53},
  {0xc04d893d,"-P 2ab 2"        ,"P m a n"   ,' ', 53},
  {0xb36aed9f,"-P 2a 2ac"       ,"P c c a"   ,' ', 54},
  {0x88dc8d86,"-P 2b 2c"        ,"P c c b"   ,' ', 54},
  {0x25d8ca19,"-P 2a 2b"        ,"P b a a"   ,' ', 54},
  {0x5bb6eda6,"-P 2ac 2c"       ,"P c a a"   ,' ', 54},
  {0xd3456635,"-P 2bc 2b"       ,"P b c b"   ,' ', 54},
  {0x6f49ae97,"-P 2b 2ab"       ,"P b a b"   ,' ', 54},
  {0xa3851163,"-P 2 2ab"        ,"P b a m"   ,' ', 55},
  {0x8b3b9e72,"-P 2bc 2"        ,"P m c b"   ,' ', 55},
  {0x364e3ba9,"-P 2ac 2ac"      ,"P c m a"   ,' ', 55},
  {0x0e8156fc,"-P 2ab 2ac"      ,"P c c n"   ,' ', 56},
  {0x31173243,"-P 2ac 2bc"      ,"P n a a"   ,' ', 56},
  {0xbebdb03a,"-P 2bc 2ab"      ,"P b n b"   ,' ', 56},
  {0x3a7e15cd,"-P 2c 2b"        ,"P b c m"   ,' ', 57},
  {0xddeb36dc,"-P 2c 2ac"       ,"P c a m"   ,' ', 57},
  {0xe45d36f0,"-P 2ac 2a"       ,"P m c a"   ,' ', 57},
  {0x46105247,"-P 2b 2a"        ,"P m a b"   ,' ', 57},
  {0x48201c16,"-P 2a 2ab"       ,"P b m a"   ,' ', 57},
  {0x280f97bc,"-P 2bc 2c"       ,"P c m b"   ,' ', 57},
  {0x43493b98,"-P 2 2n"         ,"P n n m"   ,' ', 58},
  {0x609e9307,"-P 2n 2"         ,"P m n n"   ,' ', 58},
  {0xc4f39323,"-P 2n 2n"        ,"P n m n"   ,' ', 58},
  {0x57337891,"P 2 2ab -1ab"    ,"P m m n"   ,'1', 59},
  {0xadb55f32,"-P 2ab 2a"       ,"P m m n"   ,'2', 59},
  {0x2282419e,"P 2bc 2 -1bc"    ,"P n m m"   ,'1', 59},
  {0xdab23f36,"-P 2c 2bc"       ,"P n m m"   ,'2', 59},
  {0xe9d1ae0a,"P 2ac 2ac -1ac"  ,"P m n m"   ,'1', 59},
  {0x0ff83b85,"-P 2c 2a"        ,"P m n m"   ,'2', 59},
  {0x5518bd4f,"-P 2n 2ab"       ,"P b c n"   ,' ', 60},
  {0xc3aa9ac9,"-P 2n 2c"        ,"P c a n"   ,' ', 60},
  {0xa8ec36ed,"-P 2a 2n"        ,"P n c a"   ,' ', 60},
  {0x2f569e56,"-P 2bc 2n"       ,"P n a b"   ,' ', 60},
  {0xd1db18b8,"-P 2ac 2b"       ,"P b n a"   ,' ', 60},
  {0xe5245b89,"-P 2b 2ac"       ,"P c n b"   ,' ', 60},
  {0xbc23ceb7,"-P 2ac 2ab"      ,"P b c a"   ,' ', 61},
  {0x45f741b3,"-P 2bc 2ac"      ,"P c a b"   ,' ', 61},
  {0x5cefe44c,"-P 2ac 2n"       ,"P n m a"   ,' ', 62},
  {0xe6c3487d,"-P 2bc 2a"       ,"P m n b"   ,' ', 62},
  {0x5786c3c2,"-P 2c 2ab"       ,"P b n m"   ,' ', 62},
  {0xae524cc6,"-P 2n 2ac"       ,"P c m n"   ,' ', 62},
  {0x0d664508,"-P 2n 2a"        ,"P m c n"   ,' ', 62},
  {0xb74ae939,"-P 2c 2n"        ,"P n a m"   ,' ', 62},
  {0xc0dd6737,"-C 2c 2"         ,"C m c m"   ,' ', 63},
  {0xa4e78801,"-C 2c 2c"        ,"C c m m"   ,' ', 63},
  {0x637ac153,"-A 2a 2a"        ,"A m m a"   ,' ', 63},
  {0x7af93f9c,"-A 2 2a"         ,"A m a m"   ,' ', 63},
  {0x5ea5920d,"-B 2 2b"         ,"B b m m"   ,' ', 63},
  {0xdf056b84,"-B 2b 2"         ,"B m m b"   ,' ', 63},
  {0xd95e99f8,"-C 2bc 2"        ,"C m c a"   ,' ', 64},
  {0x014f7f4e,"-C 2bc 2bc"      ,"C c m b"   ,' ', 64},
  {0xe2da38da,"-A 2ac 2ac"      ,"A b m a"   ,' ', 64},
  {0x3e8a450b,"-A 2 2ac"        ,"A c a m"   ,' ', 64},
  {0xe28e9b8d,"-B 2 2bc"        ,"B b c m"   ,' ', 64},
  {0xc686954b,"-B 2bc 2"        ,"B m a b"   ,' ', 64},
  {0xcc263714,"-C 2 2"          ,"C m m m"   ,' ', 65},
  {0xc6d2361c,"-A 2 2"          ,"A m m m"   ,' ', 65},
  {0x1ad6e89a,"-B 2 2"          ,"B m m m"   ,' ', 65},
  {0xa81cd822,"-C 2 2c"         ,"C c c m"   ,' ', 66},
  {0xdf51c8d3,"-A 2a 2"         ,"A m a a"   ,' ', 66},
  {0x9b761113,"-B 2b 2b"        ,"B b m b"   ,' ', 66},
  {0xd5a5c9db,"-C 2b 2"         ,"C m m a"   ,' ', 67},
  {0x698ec05b,"-C 2b 2b"        ,"C m m b"   ,' ', 67},
  {0x4772cf95,"-A 2c 2c"        ,"A b m m"   ,' ', 67},
  {0x82a14c8b,"-A 2 2c"         ,"A c m m"   ,' ', 67},
  {0xa6fde11a,"-B 2 2c"         ,"B m c m"   ,' ', 67},
  {0x03551655,"-B 2c 2"         ,"B m a m"   ,' ', 67},
  {0xd12ab103,"C 2 2 -1bc"      ,"C c c a"   ,'1', 68},
  {0x0db42f6d,"-C 2b 2bc"       ,"C c c a"   ,'2', 68},
  {0xd12ab103,"C 2 2 -1bc"      ,"C c c b"   ,'1', 68},
  {0xb19f26ed,"-C 2b 2c"        ,"C c c b"   ,'2', 68},
  {0x6a64cf1d,"A 2 2 -1ac"      ,"A b a a"   ,'1', 68},
  {0x9b22b244,"-A 2a 2c"        ,"A b a a"   ,'2', 68},
  {0x6a64cf1d,"A 2 2 -1ac"      ,"A c a a"   ,'1', 68},
  {0x5ef1315a,"-A 2ac 2c"       ,"A c a a"   ,'2', 68},
  {0xb660119b,"B 2 2 -1bc"      ,"B b c b"   ,'1', 68},
  {0x82f5efdc,"-B 2bc 2b"       ,"B b c b"   ,'2', 68},
  {0xb660119b,"B 2 2 -1bc"      ,"B b a b"   ,'1', 68},
  {0x275d1893,"-B 2b 2bc"       ,"B b a b"   ,'2', 68},
  {0xcd91580e,"-F 2 2"          ,"F m m m"   ,' ', 69},
  {0xd31c4512,"F 2 2 -1d"       ,"F d d d"   ,'1', 70},
  {0xc3fdc94f,"-F 2uv 2vw"      ,"F d d d"   ,'2', 70},
  {0x4a8ca5fc,"-I 2 2"          ,"I m m m"   ,' ', 71},
  {0xb2d4d6eb,"-I 2 2c"         ,"I b a m"   ,' ', 72},
  {0x530f5b33,"-I 2a 2"         ,"I m c b"   ,' ', 72},
  {0xcb2c5c75,"-I 2b 2b"        ,"I c m a"   ,' ', 72},
  {0x770755f5,"-I 2b 2c"        ,"I b c a"   ,' ', 73},
  {0x177c21a4,"-I 2a 2b"        ,"I c a b"   ,' ', 73},
  {0x8f5f26e2,"-I 2b 2"         ,"I m m a"   ,' ', 74},
  {0xef2452b3,"-I 2a 2a"        ,"I m m b"   ,' ', 74},
  {0x6e84ab3a,"-I 2c 2c"        ,"I b m m"   ,' ', 74},
  {0x0effdf6b,"-I 2 2b"         ,"I c m m"   ,' ', 74},
  {0xf6a7ac7c,"-I 2 2a"         ,"I m c m"   ,' ', 74},
  {0x96dcd82d,"-I 2c 2"         ,"I m a m"   ,' ', 74},
  {0xe194247e,"P 4"             ,"P 4"       ,' ', 75},
  {0xeb09f6d0,"P 4w"            ,"P 41"      ,' ', 76},
  {0x1e1ef133,"P 4c"            ,"P 42"      ,' ', 77},
  {0x21968a94,"P 4cw"           ,"P 43"      ,' ', 78},
  {0x8d2db265,"I 4"             ,"I 4"       ,' ', 79},
  {0x04c2f7b9,"I 4bw"           ,"I 41"      ,' ', 80},
  {0x5b1380ec,"P -4"            ,"P -4"      ,' ', 81},
  {0x429018fc,"I -4"            ,"I -4"      ,' ', 82},
  {0x7f473453,"-P 4"            ,"P 4/m"     ,' ', 83},
  {0x07df08e7,"-P 4c"           ,"P 42/m"    ,' ', 84},
  {0x8d9739ab,"P 4ab -1ab"      ,"P 4/n"     ,'1', 85},
  {0x28870b39,"-P 4a"           ,"P 4/n"     ,'2', 85},
  {0xda4091d2,"P 4n -1n"        ,"P 42/n"    ,'1', 86},
  {0xcfc232bb,"-P 4bc"          ,"P 42/n"    ,'2', 86},
  {0xe426791e,"-I 4"            ,"I 4/m"     ,' ', 87},
  {0xdd05fe66,"I 4bw -1bw"      ,"I 41/a"    ,'1', 88},
  {0xf5a09fa3,"-I 4ad"          ,"I 41/a"    ,'2', 88},
  {0x781566a5,"P 4 2"           ,"P 4 2 2"   ,' ', 89},
  {0x529dbafe,"P 4ab 2ab"       ,"P 4 21 2"  ,' ', 90},
  {0x1129373c,"P 4w 2c"         ,"P 41 2 2"  ,' ', 91},
  {0x7abbc4b0,"P 4abw 2nw"      ,"P 41 21 2" ,' ', 92},
  {0x00a13651,"P 4c 2"          ,"P 42 2 2"  ,' ', 93},
  {0xb11a4eb4,"P 4n 2n"         ,"P 42 21 2" ,' ', 94},
  {0xc507cd54,"P 4cw 2c"        ,"P 43 2 2"  ,' ', 95},
  {0x5af2c5d9,"P 4nw 2abw"      ,"P 43 21 2" ,' ', 96},
  {0xcc0e6dad,"I 4 2"           ,"I 4 2 2"   ,' ', 97},
  {0xac45c901,"I 4bw 2bw"       ,"I 41 2 2"  ,' ', 98},
  {0x93373957,"P 4 -2"          ,"P 4 m m"   ,' ', 99},
  {0xd7961388,"P 4 -2ab"        ,"P 4 b m"   ,' ',100},
  {0x53dd13c8,"P 4c -2c"        ,"P 42 c m"  ,' ',101},
  {0x25fb5987,"P 4n -2n"        ,"P 42 n m"  ,' ',102},
  {0x4124340e,"P 4 -2c"         ,"P 4 c c"   ,' ',103},
  {0x375a3973,"P 4 -2n"         ,"P 4 n c"   ,' ',104},
  {0x81ce1e91,"P 4c -2"         ,"P 42 m c"  ,' ',105},
  {0x5b8d5c98,"P 4c -2ab"       ,"P 42 b c"  ,' ',106},
  {0xd27b6f7f,"I 4 -2"          ,"I 4 m m"   ,' ',107},
  {0x2a231c68,"I 4 -2c"         ,"I 4 c m"   ,' ',108},
  {0x6ee4a2b9,"I 4bw -2"        ,"I 41 m d"  ,' ',109},
  {0xf0719502,"I 4bw -2c"       ,"I 41 c d"  ,' ',110},
  {0xae367c74,"P -4 2"          ,"P -4 2 m"  ,' ',111},
  {0xc47b0b46,"P -4 2c"         ,"P -4 2 c"  ,' ',112},
  {0x84bea02f,"P -4 2ab"        ,"P -4 21 m" ,' ',113},
  {0x67395465,"P -4 2n"         ,"P -4 21 c" ,' ',114},
  {0xd50d13d5,"P -4 -2"         ,"P -4 m 2"  ,' ',115},
  {0x3a48de91,"P -4 -2c"        ,"P -4 c 2"  ,' ',116},
  {0x00b4acac,"P -4 -2ab"       ,"P -4 b 2"  ,' ',117},
  {0xf7b01eef,"P -4 -2n"        ,"P -4 n 2"  ,' ',118},
  {0x7a127a20,"I -4 -2"         ,"I -4 m 2"  ,' ',119},
  {0xc66b145c,"I -4 -2c"        ,"I -4 c 2"  ,' ',120},
  {0x93e9bb9a,"I -4 2"          ,"I -4 2 m"  ,' ',121},
  {0xf3a21f36,"I -4 2bw"        ,"I -4 2 d"  ,' ',122},
  {0xb8f2113e,"-P 4 2"          ,"P 4/m m m" ,' ',123},
  {0x2de869d5,"-P 4 2c"         ,"P 4/m c c" ,' ',124},
  {0xb7fa689e,"P 4 2 -1ab"      ,"P 4/n b m" ,'1',125},
  {0x026e1f11,"-P 4a 2b"        ,"P 4/n b m" ,'2',125},
  {0x08a4eb7f,"P 4 2 -1n"       ,"P 4/n n c" ,'1',126},
  {0x51368ca3,"-P 4a 2bc"       ,"P 4/n n c" ,'2',126},
  {0xd5a84150,"-P 4 2ab"        ,"P 4/m b m" ,' ',127},
  {0xf1b4f061,"-P 4 2n"         ,"P 4/m n c" ,' ',128},
  {0x81d6f60a,"P 4ab 2ab -1ab"  ,"P 4/n m m" ,'1',129},
  {0x6f344f7f,"-P 4a 2a"        ,"P 4/n m m" ,'2',129},
  {0x9aee9967,"P 4ab 2n -1ab"   ,"P 4/n c c" ,'1',130},
  {0x8d6a1517,"-P 4a 2ac"       ,"P 4/n c c" ,'2',130},
  {0x5a9858b3,"-P 4c 2"         ,"P 42/m m c",' ',131},
  {0xcf822058,"-P 4c 2c"        ,"P 42/m c m",' ',132},
  {0x00e9903a,"P 4n 2c -1n"     ,"P 42/n b c",'1',133},
  {0xf9df795b,"-P 4ac 2b"       ,"P 42/n b c",'2',133},
  {0xaad7368d,"P 4n 2 -1n"      ,"P 42/n n m",'1',134},
  {0xaa87eae9,"-P 4ac 2bc"      ,"P 42/n n m",'2',134},
  {0xa5c2e426,"-P 4c 2ab"       ,"P 42/m b c",' ',135},
  {0x53c72d93,"-P 4n 2n"        ,"P 42/m n m",' ',136},
  {0xe391d7d2,"P 4n 2n -1n"     ,"P 42/n m c",'1',137},
  {0x0685c5ce,"-P 4ac 2a"       ,"P 42/n m c",'2',137},
  {0xf8a9b8bf,"P 4n 2ab -1n"    ,"P 42/n c m",'1',138},
  {0xe4db9fa6,"-P 4ac 2ac"      ,"P 42/n c m",'2',138},
  {0x7a639eef,"-I 4 2"          ,"I 4/m m m" ,' ',139},
  {0xd7e16f02,"-I 4 2c"         ,"I 4/m c m" ,' ',140},
  {0xe12277ea,"I 4bw 2bw -1bw"  ,"I 41/a m d",'1',141},
  {0xb96409e0,"-I 4bd 2"        ,"I 41/a m d",'2',141},
  {0x8bdd8359,"I 4bw 2aw -1bw"  ,"I 41/a c d",'1',142},
  {0x08d64e97,"-I 4bd 2c"       ,"I 41/a c d",'2',142},
  {0x329c980e,"P 3"             ,"P 3"       ,' ',143},
  {0x86dd212b,"P 31"            ,"P 31"      ,' ',144},
  {0xd70a5f46,"P 32"            ,"P 32"      ,' ',145},
  {0x8e4e25f8,"R 3"             ,"R 3"       ,'H',146},
  {0xd5a0aa2d,"P 3*"            ,"R 3"       ,'R',146},
  {0xfdd759b5,"-P 3"            ,"P -3"      ,' ',147},
  {0xbe8d0d7f,"-R 3"            ,"R -3"      ,'H',148},
  {0xd9a29bac,"-P 3*"           ,"R -3"      ,'R',148},
  {0x65b7a72b,"P 3 2"           ,"P 3 1 2"   ,' ',149},
  {0xc1840a7a,"P 3 2\""         ,"P 3 2 1"   ,' ',150},
  {0x97e2dfd5,"P 31 2c (0 0 1)" ,"P 31 1 2"  ,' ',151},
  {0x33d17284,"P 31 2\""        ,"P 31 2 1"  ,' ',152},
  {0xe39f36b4,"P 32 2c (0 0 -1)","P 32 1 2"  ,' ',153},
  {0x47ac9be5,"P 32 2\""        ,"P 32 2 1"  ,' ',154},
  {0x46ebee09,"R 3 2\""         ,"R 3 2"     ,'H',155},
  {0xa20b8591,"P 3* 2"          ,"R 3 2"     ,'R',155},
  {0x9f4cffaa,"P 3 -2\""        ,"P 3 m 1"   ,' ',156},
  {0x39859b12,"P 3 -2"          ,"P 3 1 m"   ,' ',157},
  {0xe04fe588,"P 3 -2\"c"       ,"P 3 c 1"   ,' ',158},
  {0xec0db0dd,"P 3 -2c"         ,"P 3 1 c"   ,' ',159},
  {0x6e32558c,"R 3 -2\""        ,"R 3 m"     ,'H',160},
  {0xb951b4f7,"P 3* -2"         ,"R 3 m"     ,'R',160},
  {0xf7d1a830,"R 3 -2\"c"       ,"R 3 c"     ,'H',161},
  {0x219be015,"P 3* -2n"        ,"R 3 c"     ,'R',161},
  {0xf74c7f83,"-P 3 2"          ,"P -3 1 m"  ,' ',162},
  {0x69dc2f41,"-P 3 2c"         ,"P -3 1 c"  ,' ',163},
  {0xfc3edafb,"-P 3 2\""        ,"P -3 m 1"  ,' ',164},
  {0xe9d82a99,"-P 3 2\"c"       ,"P -3 c 1"  ,' ',165},
  {0x6df507a9,"-R 3 2\""        ,"R -3 m"    ,'H',166},
  {0x1c80e47a,"-P 3* 2"         ,"R -3 m"    ,'R',166},
  {0x9a7d09d3,"-R 3 2\"c"       ,"R -3 c"    ,'H',167},
  {0xbb691c91,"-P 3* 2n"        ,"R -3 c"    ,'R',167},
  {0xa2ddaf47,"P 6"             ,"P 6"       ,' ',168},
  {0x81a2968a,"P 61"            ,"P 61"      ,' ',169},
  {0x84b83b8b,"P 65"            ,"P 65"      ,' ',170},
  {0x62d97e5d,"P 62"            ,"P 62"      ,' ',171},
  {0xb22fefa7,"P 64"            ,"P 64"      ,' ',172},
  {0xdddeb565,"P 6c"            ,"P 63"      ,' ',173},
  {0x35eb6ccb,"P -6"            ,"P -6"      ,' ',174},
  {0x32dacfb6,"-P 6"            ,"P 6/m"     ,' ',175},
  {0x06c1ae99,"-P 6c"           ,"P 63/m"    ,' ',176},
  {0x90230e5f,"P 6 2"           ,"P 6 2 2"   ,' ',177},
  {0xeed642a2,"P 61 2 (0 0 -1)" ,"P 61 2 2"  ,' ',178},
  {0x891aeb21,"P 65 2 (0 0 1)"  ,"P 65 2 2"  ,' ',179},
  {0x6f55448f,"P 62 2c (0 0 1)" ,"P 62 2 2"  ,' ',180},
  {0x655e4cd9,"P 64 2c (0 0 -1)","P 64 2 2"  ,' ',181},
  {0xa4386f70,"P 6c 2c"         ,"P 63 2 2"  ,' ',182},
  {0x3b0a2d17,"P 6 -2"          ,"P 6 m m"   ,' ',183},
  {0xb74f3653,"P 6 -2c"         ,"P 6 c c"   ,' ',184},
  {0xf4d8e94d,"P 6c -2"         ,"P 63 c m"  ,' ',185},
  {0x789df209,"P 6c -2c"        ,"P 63 m c"  ,' ',186},
  {0xeeacb736,"P -6 2"          ,"P -6 m 2"  ,' ',187},
  {0xfb4a4754,"P -6c 2"         ,"P -6 c 2"  ,' ',188},
  {0x55b5be6a,"P -6 -2"         ,"P -6 2 m"  ,' ',189},
  {0xcb25eea8,"P -6c -2c"       ,"P -6 2 c"  ,' ',190},
  {0xf1fc7952,"-P 6 2"          ,"P 6/m m m" ,' ',191},
  {0x87f9e8ca,"-P 6 2c"         ,"P 6/m c c" ,' ',192},
  {0x1eb41bd9,"-P 6c 2"         ,"P 63/m c m",' ',193},
  {0x68b18a41,"-P 6c 2c"        ,"P 63/m m c",' ',194},
  {0x5843870d,"P 2 2 3"         ,"P 2 3"     ,' ',195},
  {0x93e38c71,"F 2 2 3"         ,"F 2 3"     ,' ',196},
  {0xdc4003c1,"I 2 2 3"         ,"I 2 3"     ,' ',197},
  {0xf9e6a645,"P 2ac 2ab 3"     ,"P 21 3"    ,' ',198},
  {0x7ec4457b,"I 2b 2c 3"       ,"I 21 3"    ,' ',199},
  {0x72e55913,"-P 2 2 3"        ,"P m -3"    ,' ',200},
  {0x265ce726,"P 2 2 3 -1n"     ,"P n -3"    ,'1',201},
  {0xd419d8a3,"-P 2ab 2bc 3"    ,"P n -3"    ,'2',201},
  {0x58320b8d,"-F 2 2 3"        ,"F m -3"    ,' ',202},
  {0x7de7c89b,"F 2 2 3 -1d"     ,"F d -3"    ,'1',203},
  {0x159ca8d3,"-F 2uv 2vw 3"    ,"F d -3"    ,'2',203},
  {0xe23893df,"-I 2 2 3"        ,"I m -3"    ,' ',204},
  {0x1d5f4d3f,"-P 2ac 2ab 3"    ,"P a -3"    ,' ',205},
  {0x7ce42b66,"-I 2b 2c 3"      ,"I a -3"    ,' ',206},
  {0x93a6edeb,"P 4 2 3"         ,"P 4 3 2"   ,' ',207},
  {0xc55ae72a,"P 4n 2 3"        ,"P 42 3 2"  ,' ',208},
  {0x25d65cf1,"F 4 2 3"         ,"F 4 3 2"   ,' ',209},
  {0xddc15ef3,"F 4d 2 3"        ,"F 41 3 2"  ,' ',210},
  {0x014b7ed2,"I 4 2 3"         ,"I 4 3 2"   ,' ',211},
  {0x4bca2b34,"P 4acd 2ab 3"    ,"P 43 3 2"  ,' ',212},
  {0x0e761d2d,"P 4bd 2ab 3"     ,"P 41 3 2"  ,' ',213},
  {0xbb92b652,"I 4bd 2c 3"      ,"I 41 3 2"  ,' ',214},
  {0x6fc31c7b,"P -4 2 3"        ,"P -4 3 m"  ,' ',215},
  {0x09abdebe,"F -4 2 3"        ,"F -4 3 m"  ,' ',216},
  {0x75d77c69,"I -4 2 3"        ,"I -4 3 m"  ,' ',217},
  {0x393f16ba,"P -4n 2 3"       ,"P -4 3 n"  ,' ',218},
  {0x0968710f,"F -4c 2 3"       ,"F -4 3 c"  ,' ',219},
  {0xe2ed982e,"I -4bd 2c 3"     ,"I -4 3 d"  ,' ',220},
  {0x74c407d3,"-P 4 2 3"        ,"P m -3 m"  ,' ',221},
  {0xc5e3dd9f,"P 4 2 3 -1n"     ,"P n -3 n"  ,'1',222},
  {0xc7e69ab6,"-P 4a 2bc 3"     ,"P n -3 n"  ,'2',222},
  {0xbabd71c4,"-P 4n 2 3"       ,"P m -3 n"  ,' ',223},
  {0x0b9aab88,"P 4n 2 3 -1n"    ,"P n -3 m"  ,'1',224},
  {0xdbdb6bfc,"-P 4bc 2bc 3"    ,"P n -3 m"  ,'2',224},
  {0x5225b614,"-F 4 2 3"        ,"F m -3 m"  ,' ',225},
  {0x481e9f10,"-F 4c 2 3"       ,"F m -3 c"  ,' ',226},
  {0x023a8184,"F 4d 2 3 -1d"    ,"F d -3 m"  ,'1',227},
  {0x38f5e458,"-F 4vw 2vw 3"    ,"F d -3 m"  ,'2',227},
  {0x82b6e512,"F 4d 2 3 -1cd"   ,"F d -3 c"  ,'1',228},
  {0x4d325ad3,"-F 4cvw 2vw 3"   ,"F d -3 c"  ,'2',228},
  {0x52045e84,"-I 4 2 3"        ,"I m -3 m"  ,' ',229},
  {0x407de7c1,"-I 4bd 2c 3"     ,"I a -3 d"  ,' ',230}
};
int sgdata_size = sizeof( sgdata ) / sizeof( sgdata[0] );


} // namespace data


} // namespace clipper
