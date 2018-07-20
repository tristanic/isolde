//  $Id: mmdb_math_bfgsmin.h $
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
//  **** Module  :  BFGSMin  <interface>
//       ~~~~~~~~~
//  **** Classes :  mmdb::math::BFGSMin  ( minimization driver )
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#ifndef  __MMDB_MATH_BFGSMin__
#define  __MMDB_MATH_BFGSMin__

#include <stdlib.h>

#include "mmdb_mattype.h"
#include "mmdb_math_linalg.h"


namespace mmdb  {

  namespace math  {

    //  ==============================================================

    enum BFGS_RC  {
      BFGS_TooFewDigits       = -2,
      BFGS_WrongSpaceDim      = -1,
      BFGS_NoTermination      =  0,
      BFGS_SmallGradient      =  1,
      BFGS_SmallStep          =  2,
      BFGS_LineSearchComplete =  3,
      BFGS_IterationLimit     =  4,
      BFGS_LargeSteps         =  5,
      BFGS_Stopped            =  6
    };

    typedef void BFGSMinFunc ( void * UserData, int N, rvector X,
                               realtype & F );
    typedef BFGSMinFunc * PBFGSMinFunc;

    typedef void BFGSPrintFunc ( void * UserData, int N, int Itn,
                                 rvector X, rvector G, realtype F );
    typedef BFGSPrintFunc * PBFGSPrintFunc;

    DefineClass(BFGSMin);

    class BFGSMin  {

      public :

        BFGSMin ();
        virtual ~BFGSMin();

        virtual void  MinFunc  ( rvector X, realtype & F );
        virtual void  Print    ( int Itn, rvector X, rvector G,
                                 realtype F );

        void  SetMinFunction   ( void * UserData, PBFGSMinFunc   Fnc );
        void  SetPrintFunction ( void * UserData, PBFGSPrintFunc Fnc );


        // ======================================================
        //
        //    .--------------------------------------------.
        //    |                                            |
        //    |     UNCONSTRAINED MINIMIZATION DRIVER      |
        //    |                                            |
        //    `--------------------------------------------'
        //
        //    Finds a minimum of function F(X), X is vector [1..N],
        //  defined by virtual MinFunc. Virtual Print provides
        //  information on every iteration step.
        //
        //
        //               Input  parameters  :
        //             -----------------------
        //
        //    N        is the dimension the minimization space
        //
        //    x0       [1..N] is the initial point for minimization
        //
        //    TypX     [1..N] is the array of the typical ranges of
        //          X - components,  which are used for the scaling.
        //          If  TypX<=0.0  then  1.0  will be substituted
        //
        //    Digits   is the number of valid decimal digits in
        //          the calculated value of minimizing function ( F ).
        //          If  Digits<=0  then the Driver will consider
        //          that the  F  is computed with usual machine's
        //          noise
        //
        //    ItnLmt   is the maximum available number of iterations.
        //          If  ItnLmt=0  then  100  will be substituted
        //
        //    TypF     is the expected absolute value of  F  in the
        //          minimum,  which is used in the stop criterion.
        //          If  TypF<=0.0  then  1.0  will be substituted
        //
        //    GrdTol   is the desired absolute value of the gradient
        //          vector in the minimum of  F .  If  GrdTol<=0.0
        //          then the some value correlated with machine's
        //          noise will be substituted
        //
        //    StpTol   is the minimum available step for the minimi-
        //          zation.  The execution stops if the distance
        //          between two consequential approximation will be
        //          less then  StpTol .  If  StpTol<=0.0  then the
        //          some value correlated with machine's  noise
        //          will be substituted
        //
        //    MaxStp   is the maximum available step for then minimi-
        //          zation.  This parameter only prevents the appea-
        //          rance of the too large steps,  but the execution
        //          stops if more than  5  steps with length of MaxStep
        //          will consequently appear.
        //
        //
        //
        //               Outpute  parameters  :
        //             --------------------------
        //
        //    x0        will be the point at which the minimisation
        //          had stopped
        //
        //    Func      will be the function's value at  x0
        //
        //    TermCode  will be the reason of stopping :
        //
        //         1 <=>  the norm of gradient vector at  x0  is
        //               less than  GrdTol ;  the  x0  is probable
        //               point of the minimum
        //         2 <=>  the distance between two last approxima-
        //               tions was less than  StpTol ;  the  x0
        //               may be the point of minimum
        //         3 <=>  the gradient length is greater than
        //               GrdTol ,  but future minimization fails ;
        //               it may be the consequence of the errors
        //               at the computing of gradient, but also
        //               x0 could be the point of minimum
        //         4 <=>  the iteration limit had been exchausted
        //         5 <=>  more than  5  steps with length of
        //               MaxStp  had been made
        //         6 <=>  the termination key ( Esc or End )
        //               had been pressed.
        //
        //
        // ========================================================

        void  BFGS_Driver      ( int        MinN,
                                 rvector    x0,
                                 rvector    TypX,
                                 realtype & FuncValue,
                                 int      & TerminationCode,
                                 int        Digits   = 0,
                                 int        ItnLmt   = 0,
                                 realtype   TypF     = 0.0,
                                 realtype   GrdTol   = 0.0,
                                 realtype   StpTol   = 0.0,
                                 realtype   MaxStp   = MaxReal,
                                 bool       Hess     = false,
                                 rvector    LowLimit = NULL,
                                 rvector    TopLimit = NULL );

        void  Stop();  // generates stop signal to stop optimization


      protected :

        PBFGSMinFunc    MFunc;
        void *          MFuncData;
        PBFGSPrintFunc  PFunc;
        void *          PFuncData;

        int             N,NAlloc;
        rmatrix         Hsn;
        rvector         TL,LL,XOpt,XPlus,Sx,SN,HDiag,GradX,GPlus;
        rvector         StepSize,FNeighbor;
        rvector         us,uy,ut;
        bvector         Freese;
        realtype        Func,FPlus,FOpt;
        realtype        TakenLambda;
        bool            ForDiff;  // if True then forward differences are
                                  // used for the 1st estimation of the
                                  // Hessian (which is less expensive),
                                  // otherwise central differences will
                                  // be employed (which is more expensive).
        bool            CalcHess;

        realtype        Etha,SqrtEtha,CubertEtha,TpF,GrdEps,StpEps,MxStep;
        realtype        SqrtEps;
        int             CnsMax,MaxItn,TermCode;
        bool            ModF;

        void  MinFunc1      ( rvector X, realtype & F );
        void  UMInCk        ( rvector  x0,     rvector  TypX,
                              int      Digits, realtype TypF,
                              realtype GrdTol, realtype StpTol,
                              realtype MaxStp, int      ItnLmt );
        void  UMStop0       ( rvector x0, rvector Grad );
        void  UMStop        ( rvector x0, rvector Grad, int RetCode,
                              int ItnCnt, bool MaxTkn );

        virtual void Gradient ( rvector X, rvector G, realtype Fc );
        virtual void FDHessF  ( realtype Fc, rvector X );

        void  FDGrad        ( rvector X, rvector G, realtype Fc );
        void  CDGrad        ( rvector X, rvector G );
        void  MdHess        ( rmatrix H, rvector HDg );
        void  InitHessUnFac ( realtype F,  rmatrix H );
        void  BFGSUnFac     ( rvector  Xc,      rvector Xp,
                              rvector  Gc,      rvector Gp,
                              bool     AnalGrad, rvector HDg,
                              rmatrix  H );
        void  Choose_Lambda ( rvector X, rvector S, realtype & Lambda0 );
        void  LineSearch    ( rvector    px0,       rvector   G,
                              rvector    P,         realtype pFunc,
                              int     & RetCode,    bool & MaxTkn );
        void  GetMemory     ();
        void  FreeMemory    ();
        void  Relax         ();
        void  CopyPlus      ( rvector x0 );

    };

  }

}

#endif


