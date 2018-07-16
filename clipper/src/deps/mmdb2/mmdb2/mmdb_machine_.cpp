//  $Id: mmdb_machine.cpp $
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
//  **** Module  :   Machine  <interface>
//       ~~~~~~~~~
//  **** Functions : mmdb::machine::GetMachineID   - returns ID code
//       ~~~~~~~~~~~                                 for the machine
//                   mmdb::machine::GetMachineName - returns name of
//                                                   the machine
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#include "mmdb_machine_.h"

namespace mmdb  {

  namespace machine  {

#ifdef CALL_LIKE_SUN

    int GetMachineID()  {
    int k = CALL_LIKE_SUN;
      switch (k)  {
        case 1  : return MACHINE_ALLIANT;
        case 2  : return MACHINE_CONVEX;
        case 3  : return MACHINE_ESV;
        case 4  : return MACHINE_SGI;
        case 5  : return MACHINE_SOLBOURNE;
        case 6  : return MACHINE_SOLARIS;
        case 7  : return MACHINE_ALPHA;
        case 8  : return MACHINE_F2C_G77;
        case 9  : return MACHINE_LINUX;
        default : return MACHINE_UNKNOWN;
      }
    }

#elif defined(CALL_LIKE_HPUX)

    int GetMachineID()  {
    int k = CALL_LIKE_HPUX;
      switch (k)  {
        case 1  : return MACHINE_RS6000;
        case 2  : return MACHINE_HP9000;
        default : return MACHINE_UNKNOWN;
      }
    }

#elif defined(CALL_LIKE_STARDENT)

    int GetMachineID()  {
    int k = CALL_LIKE_STARDENT;
      switch (k)  {
        case 1  : return MACHINE_ARDENT;
        case 2  : return MACHINE_TITAN;
        case 3  : return MACHINE_STARDENT;
        default : return MACHINE_UNKNOWN;
      }
    }

#elif defined(CALL_LIKE_VMS)

    int GetMachineID()  {
      return MACHINE_VMS;
    }

#elif defined(CALL_LIKE_MVS)

    int GetMachineID()  {
      return MACHINE_MVS;
    }

#else

    int GetMachineID()  {
      return MACHINE_UNKNOWN;
    }

#endif

    static cpstr MCH_SGI       = cpstr("Silicon Graphics");
    static cpstr MCH_RS6000    = cpstr("IBM RS/6000");
    static cpstr MCH_ALLIANT   = cpstr("Alliant");
    static cpstr MCH_ARDENT    = cpstr("Ardent");
    static cpstr MCH_TITAN     = cpstr("Titan");
    static cpstr MCH_STARDENT  = cpstr("Stardent");
    static cpstr MCH_CONVEX    = cpstr("Convex");
    static cpstr MCH_ESV       = cpstr("Evans or Sutherland");
    static cpstr MCH_HP9000    = cpstr("Hewlett Packard 9000");
    static cpstr MCH_SOLBOURNE = cpstr("Solbourne");
    static cpstr MCH_SOLARIS   = cpstr("Solaris");
    static cpstr MCH_ALPHA     = cpstr("DEC Alpha");
    static cpstr MCH_VMS       = cpstr("A VMS machine");
    static cpstr MCH_MVS       = cpstr("MS Windows");
    static cpstr MCH_F2C_G77   = cpstr("SUN compatible");
    static cpstr MCH_LINUX     = cpstr("Linux");

    cpstr GetMachineName ( int MachineID )  {
      switch (MachineID)  {
        case MACHINE_SGI       : return MCH_SGI;
        case MACHINE_RS6000    : return MCH_RS6000;
        case MACHINE_ALLIANT   : return MCH_ALLIANT;
        case MACHINE_ARDENT    : return MCH_ARDENT;
        case MACHINE_TITAN     : return MCH_TITAN;
        case MACHINE_STARDENT  : return MCH_STARDENT;
        case MACHINE_CONVEX    : return MCH_CONVEX;
        case MACHINE_ESV       : return MCH_ESV;
        case MACHINE_HP9000    : return MCH_HP9000;
        case MACHINE_SOLBOURNE : return MCH_SOLBOURNE;
        case MACHINE_SOLARIS   : return MCH_SOLARIS;
        case MACHINE_ALPHA     : return MCH_ALPHA;
        case MACHINE_VMS       : return MCH_VMS;
        case MACHINE_MVS       : return MCH_MVS;
        case MACHINE_F2C_G77   : return MCH_F2C_G77;
        case MACHINE_LINUX     : return MCH_LINUX;
        default                :
        case MACHINE_UNKNOWN   : return pstr("Unidentified machine");
      }
    }

    cpstr GetMachineName()  {
      return GetMachineName ( GetMachineID() );
    }


  }

}
