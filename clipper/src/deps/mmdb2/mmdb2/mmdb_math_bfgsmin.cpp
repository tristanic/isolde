//  $Id: mmdb_math_bfgsmin.cpp $
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
//  **** Module  :  BFGSMin  <implementation>
//       ~~~~~~~~~
//  **** Classes :  mmdb::math::BFGSMin  ( minimization driver )
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#include <math.h>

#include "mmdb_math_bfgsmin.h"

namespace mmdb  {

  namespace math  {

    //  ==============================================================

    BFGSMin::BFGSMin()  {

      MFunc       = NULL;
      MFuncData   = NULL;
      PFunc       = NULL;
      PFuncData   = NULL;

      N           = 0;
      NAlloc      = 0;

      Hsn         = NULL;
      TL          = NULL;
      LL          = NULL;
      XOpt        = NULL;
      XPlus       = NULL;
      Sx          = NULL;
      SN          = NULL;
      HDiag       = NULL;
      GradX       = NULL;
      GPlus       = NULL;
      StepSize    = NULL;
      FNeighbor   = NULL;
      us          = NULL;
      uy          = NULL;
      ut          = NULL;
      Freese      = NULL;
      Func        = 0.0;
      FPlus       = 0.0;
      FOpt        = 0.0;
      TakenLambda = 0.0;
      ForDiff     = false;
      CalcHess    = false;

      Etha        = 0.0;
      SqrtEtha    = 0.0;
      CubertEtha  = 0.0;
      TpF         = 1.0;
      GrdEps      = 0.0;
      StpEps      = 0.0;
      MxStep      = MaxReal;
      CnsMax      = 0;
      MaxItn      = 100;
      TermCode    = BFGS_NoTermination;
      ModF        = false;

    }

    BFGSMin::~BFGSMin()  {
      FreeMemory();
    }

    void  BFGSMin::MinFunc ( rvector X, realtype & F )  {
      if (MFunc)  (*MFunc)(MFuncData,N,X,F);
            else  F = 0.0;
    }

    void  BFGSMin::MinFunc1 ( rvector X, realtype & F )  {
    int i;
      MinFunc ( X,F );
      if (ModF && (F<FOpt))  {
        for (i=1;i<=N;i++)
          XOpt[i] = X[i];
        FOpt = F;
      }
    }

    void  BFGSMin::Print ( int Itn, rvector X, rvector G, realtype F )  {
      if (PFunc)  (*PFunc)(PFuncData,N,Itn,X,G,F);
    }

    void  BFGSMin::SetMinFunction ( void * UserData, PBFGSMinFunc Fnc )  {
      MFuncData = UserData;
      MFunc     = Fnc;
    }

    void  BFGSMin::SetPrintFunction ( void * UserData, PBFGSPrintFunc Fnc )  {
      PFuncData = UserData;
      PFunc     = Fnc;
    }


    //  -------------------------------------------------------------------

    void  BFGSMin::UMInCk ( rvector  x0,      rvector  TypX,
                            int      Digits,  realtype TypF,
                            realtype GrdTol,  realtype StpTol,
                            realtype MaxStp,  int      ItnLmt )  {
    int      i;
    realtype S0,S1,S2;

      SqrtEps = sqrt(MachEps);

      if (N<1)  {
        TermCode = BFGS_WrongSpaceDim;
        return;
      }

      for (i=1;i<=N;i++)
        if (fabs(TypX[i])!=0.0)  Sx[i] = 1.0/fabs(TypX[i]);
                           else  Sx[i] = 1.0;

      if (Digits<=0)  Etha = MachEps;
      else  {
        Etha = Exp((-Digits)*log(10.0));
        if (MachEps>Etha)  Etha = MachEps;
      }
      SqrtEtha   = sqrt(Etha);
      CubertEtha = Exp ( log(Etha)/3.0 );

      if (Etha>0.01)  {
        TermCode = BFGS_TooFewDigits;
        return;
      }

      if (TypF<=0.0)  TpF = 1.0;
                else  TpF = TypF;

      S1 = Exp(log(MachEps)/3.0);
      if (GrdTol>0.0)  GrdEps = GrdTol;
      else  {
        GrdEps = sqrt(Etha);
        if (S1>GrdEps)  GrdEps = S1;
      }

      if (StpTol>0.0)  StpEps = StpTol;
                 else  StpEps = Exp ( log(MachEps)*2.0/3.0 );

      if (MaxStp>0.0)  MxStep = MaxStp;
      else  {
        S1 = 0.0;
        S2 = 0.0;
        for (i=1;i<=N;i++)  {
          S0  = Sx[i];
          S0 *= Sx[i];
          S2 += S0;
          S0 *= x0[i];
          S1 += S0*x0[i];
        }
        S1 = sqrt(S1);
        S2 = sqrt(S2);
        if (S2>S1)  MxStep = S2;
              else  MxStep = S1;
        MxStep *= 1000.0;
      }

      if (ItnLmt>0)  MaxItn = ItnLmt;
               else  MaxItn = 100;

      TermCode = BFGS_NoTermination;

    }


    //  -------------------------------------------------------------------

    void  BFGSMin::UMStop0 ( rvector x0, rvector Grad )  {
    int      i;
    realtype S,Fmax,St;

      CnsMax = 0;
      if (TpF>fabs(Func))  Fmax = TpF;
                     else  Fmax = fabs(Func);
      S = 0.0;
      for (i=1;i<=N;i++)  {
        St = fabs(x0[i]);
        if (1.0/Sx[i]>St)  St = 1.0/Sx[i];
        St = fabs(Grad[i])*St/Fmax;
        if (St>S)  S = St;
      }
      if (S>=0.001*GrdEps)  TermCode = BFGS_NoTermination;
                      else  TermCode = BFGS_SmallGradient;

    }


    //  -------------------------------------------------------------------

    void  BFGSMin::UMStop ( rvector  x0, rvector Grad,
                            int RetCode, int  ItnCnt,
                            bool MaxTkn )  {

    //  A7.2.1   :  Checking the Stop Conditions

    int      i;
    realtype Max1,Max2,MaxGrad,MaxStep, BB1,BB2;

      TermCode = BFGS_NoTermination;
      if (RetCode==1)  TermCode = BFGS_LineSearchComplete;
      else  {
        if (fabs(FPlus)>TpF)  Max2 = fabs(FPlus);
                        else  Max2 = TpF;
        MaxGrad = 0.0;
        MaxStep = 0.0;
        for (i=1;i<=N;i++)  {
          BB1 = fabs(XPlus[i]);
          BB2 = 1.0/Sx[i];
          if (BB1>BB2)  Max1 = BB1;
                  else  Max1 = BB2;
          BB1 = fabs(Grad[i])*Max1/Max2;
          if (BB1>MaxGrad) MaxGrad = BB1;
          BB2 = fabs(XPlus[i]-x0[i])/Max1;
          if (BB2>MaxStep)  MaxStep = BB2;
        }
        if      (MaxGrad<GrdEps)  TermCode = BFGS_SmallGradient;
        else if (MaxStep<StpEps)  TermCode = BFGS_SmallStep;
        else if (ItnCnt>MaxItn)   TermCode = BFGS_IterationLimit;
        else if (MaxTkn)  {
          CnsMax++;
          if (CnsMax==5)  TermCode = BFGS_LargeSteps;
        } else
          CnsMax = 0;
      }

    }


    //  -------------------------------------------------------------------

    void  BFGSMin::MdHess ( rmatrix H, rvector HDg )  {

    //  A5.5.1   :  Setting up the hessian of model

    int       i,j;
    realtype  MaxDiag,MaxOff, MinEv,Mue,MaxPosDiag;
    realtype  MaxOffl,MinDiag,MaxEv,MaxAdd,Sdd,OffRow;
    realtype  BB;

      //  Scaling
      for (i=1;i<=N;i++)
        for (j=i;j<=N;j++)
          H[i][j] /= (Sx[i]*Sx[j]);

      MaxDiag = H[1][1];
      MinDiag = H[1][1];
      MaxOff  = 0.0;
      for (i=1;i<=N;i++)  {
        if (H[i][i]>MaxDiag)  MaxDiag = H[i][i];
        if (H[i][i]<MinDiag)  MinDiag = H[i][i];
        if (i<N)
          for (j=i+1;j<=N;j++)  {
            BB = fabs(H[i][j]);
            if (BB>MaxOff)  MaxOff = BB;
          }
      }
      MaxPosDiag = 0.0;
      if (MaxDiag>MaxPosDiag)  MaxPosDiag = MaxDiag;

      //  Computing the shift of the spectra (the  Mue)
      if (MinDiag>SqrtEps*MaxPosDiag)  Mue = 0.0;
      else  {
        Mue      = 2.0*(MaxPosDiag-MinDiag)*SqrtEps-MinDiag;
        MaxDiag += Mue;
      }
      BB = MaxOff*(1.0+2.0*SqrtEps);
      if (BB>MaxDiag)  {
        Mue     = Mue+(MaxOff-MaxDiag)+2.0*SqrtEps*MaxOff;
        MaxDiag = BB;
      }
      if (MaxDiag==0.0)  {  //  H = 0
        Mue     = 1.0;
        MaxDiag = 1.0;
      }
      if (Mue>0.0)
        for (i=1;i<=N;i++)
          Hsn[i][i] += Mue;

      MaxOffl = MaxOff/N;
      if (MaxDiag>MaxOffl)  MaxOffl = MaxDiag;
      MaxOffl = sqrt(MaxOffl);
      for (i=1;i<=N;i++)
        HDg[i] = H[i][i];

      PbCholDecomp ( N,HDg,MaxOffl,MachEps,H,MaxAdd );

      if (MaxAdd>0.0)  {
        MaxEv = HDg[1];
        MinEv = HDg[1];
        for (i=1;i<=N;i++)  {
          OffRow = 0.0;
          if (i>1)
            for (j=1;j<i;j++)
              OffRow += fabs(H[j][i]);
          if (i<N)
            for (j=i+1;j<=N;j++)
              OffRow += fabs(H[i][j]);
          BB = HDg[i]+OffRow;
          if (BB>MaxEv)  MaxEv = BB;
          BB = HDg[i]-OffRow;
          if (BB<MinEv)  MinEv = BB;
        }
        Sdd = (MaxEv-MinEv)*SqrtEps-MinEv;
        if (Sdd<0.0)  Sdd = 0.0;
        if (MaxAdd<Sdd)  Mue = MaxAdd;
                   else  Mue = Sdd;
        for (i=1;i<=N;i++)
          HDg[i] += Mue;

        PbCholDecomp ( N,HDg,0.0,MachEps,H,MaxAdd );

      }

      //  Scaling back
      for (i=1;i<=N;i++)  {
        if (i<N)
          for (j=i+1;j<=N;j++)
            H[i][j] *= (Sx[i]*Sx[j]);
        HDg[i] *= Sx[i]*Sx[i];
        for (j=1;j<=i;j++)
          H[i][j] *= Sx[i];
      }

    }


    //  -------------------------------------------------------------------

    void  BFGSMin::FDGrad ( rvector X, rvector G, realtype Fc )  {

    //  A5.6.4  :  Forward Finite-Differencies Approximation of
    //             the Gradient

    realtype StepSizeJ,TempJ,Fj, BB1,BB2;
    int      j;

      for (j=1;j<=N;j++)  {
        BB1 = fabs(X[j]);
        BB2 = 1.0/Sx[j];
        if (BB1>BB2)  StepSizeJ = BB1;
                else  StepSizeJ = BB2;
        if (X[j]<0.0) StepSizeJ = -StepSizeJ;
        StepSizeJ *= SqrtEtha;
        TempJ      = X[j];
        X[j]      += StepSizeJ;
        StepSizeJ  = X[j]-TempJ;
        MinFunc1 ( X,Fj );
        if (TermCode!=BFGS_NoTermination)  return;
        G[j]       = (Fj-Fc)/StepSizeJ;
        X[j]       = TempJ;
        Freese[j]  = false;
        if (TL)  {
          if ((fabs(X[j]-TL[j])<=StepSizeJ) && (G[j]<0.0))  {
            G[j] = 0.0;   Freese[j] = true;
          }
        }
        if (LL)  {
          if ((fabs(X[j]-LL[j])<=StepSizeJ) && (G[j]>0.0))  {
            G[j] = 0.0;   Freese[j] = true;
          }
        }
      }

    }


    //  -------------------------------------------------------------------

    void  BFGSMin::CDGrad ( rvector X, rvector G  )  {

    //  A5.6.4  :  Central Differencies Approximation of
    //             Gradient

    realtype  StepSizeJ,TempJ,Fp,Fm, BB1,BB2;
    int       j;

      for (j=1;j<=N;j++)  {
        BB1 = fabs(X[j]);
        BB2 = 1.0/Sx[j];
        if (BB1>BB2)  StepSizeJ = BB1;
                else  StepSizeJ = BB2;
        if (X[j]<0.0) StepSizeJ = -StepSizeJ;
        StepSizeJ *= CubertEtha;
        TempJ      = X[j];
        X[j]      += StepSizeJ;
        StepSizeJ  = X[j]-TempJ;
        MinFunc1 ( X,Fp );
        if (TermCode!=BFGS_NoTermination)  return;
        X[j]       = TempJ-StepSizeJ;
        MinFunc1 ( X,Fm );
        if (TermCode!=BFGS_NoTermination)  return;
        G[j]       = (Fp-Fm)/(2.0*StepSizeJ);
        X[j]       = TempJ;
      }

    }


    //  -------------------------------------------------------------------

    void  BFGSMin::Gradient ( rvector X, rvector G, realtype Fc )  {
      if (ForDiff)  FDGrad ( X,G,Fc );
              else  CDGrad ( X,G );
    }


    //  -------------------------------------------------------------------

    void  BFGSMin::FDHessF ( realtype Fc, rvector X )  {

    //   A5.6.2   :  Finite-Difference Approximation of
    //               the Hessian employing only the
    //               function's  values

    int       i,j;
    realtype  S,TempI,Fii,TempJ,Fij, BB1,BB2;


      for (i=1;i<=N;i++)
        if (!Freese[i])  {
          BB1 = fabs(X[i]);
          BB2 = 1.0/Sx[i];
          if (BB1>BB2)  S = BB1;
                  else  S = BB2;
          if (X[i]<0.0) S = -S;
          StepSize[i] = S*CubertEtha;
          TempI       = X[i];
          X[i]       += StepSize[i];
          StepSize[i] = X[i]-TempI;
          MinFunc1 ( X,FNeighbor[i] );
          X[i]        = TempI;
          if (TermCode!=BFGS_NoTermination)  return;
        }
      for (i=1;i<=N;i++)
        if (!Freese[i])  {
          TempI = X[i];
          X[i] += 2.0*StepSize[i];
          MinFunc1 ( X,Fii );
          if (TermCode!=BFGS_NoTermination)  return;
          Hsn[i][i] = (( Fc -FNeighbor[i] ) +
                       ( Fii-FNeighbor[i] )) /
                      (StepSize[i]*StepSize[i]);
          X[i]      = TempI+StepSize[i];
          if (i<N)
            for (j=i+1;j<=N;j++)  {
              if (!Freese[j])  {
                TempJ = X[j];
                X[j] += StepSize[j];
                MinFunc1 ( X,Fij );
                if (TermCode!=BFGS_NoTermination)  return;
                Hsn[i][j] = (( Fc -FNeighbor[i] ) +
                             ( Fij-FNeighbor[j] )) /
                            (StepSize[i]*StepSize[j]);
                X[j]      = TempJ;
              } else
                Hsn[i][j] = 0.0;
            }
          X[i] = TempI;
        } else
          for (j=i;j<=N;j++)
            Hsn[i][j] = 0.0;

    }


    //  -------------------------------------------------------------------

    void  BFGSMin::InitHessUnFac ( realtype F, rmatrix H )  {

    //  A9.4.3  -  Initialization of the unfactorized BFGS

    realtype Temp;
    int      i,j;

      Temp = fabs(F);
      if (TpF>Temp)  Temp = TpF;
      for (i=1;i<=N;i++)  {
        H[i][i] = Temp*Sx[i]*Sx[i];
        if (i<N)
          for (j=i+1;j<=N;j++)
            H[i][j] = 0.0;
      }

    }


    //  -------------------------------------------------------------------

    void  BFGSMin::BFGSUnFac (  rvector  Xc,     rvector  Xp,
                                rvector  Gc,     rvector  Gp,
                                bool  AnalGrad,  rvector  HDg,
                                rmatrix  H )  {

    //  A9.4.1  -  Calculation of the Hessian by the
    //             unfactorized BFGS

    int      i,j;
    realtype Temp1,Temp2, NormS,NormY,Tol,tt, BB;
    bool     SkipUpdate;

      Temp1 = 0.0;
      NormS = 0.0;
      NormY = 0.0;
      for (i=1;i<=N;i++)  {
        H[i][i] = HDg[i];
        us[i]   = Xp[i] - Xc[i];
        uy[i]   = Gp[i] - Gc[i];
        Temp1  += us[i]*uy[i];
        NormS  += us[i]*us[i];
        NormY  += uy[i]*uy[i];
      }

      if (Temp1>sqrt(MachEps*NormS*NormY))  {
        if (AnalGrad)  Tol = Etha;
                 else  Tol = sqrt(Etha);
        SkipUpdate = true;
        for (i=1;i<=N;i++)  {
          tt = 0.0;
          for (j=1;j<=i;j++)
            tt += H[j][i]*us[j];
          if (i<N)
            for (j=i+1;j<=N;j++)
              tt += H[i][j]*us[j];
          ut[i] = tt;
          tt    = fabs(Gc[i]);
          BB    = fabs(Gp[i]);
          if (BB>tt)  tt = BB;
          if (fabs(uy[i]-ut[i])>=Tol*tt)
            SkipUpdate = false;
        }

        if (!SkipUpdate)  {
          Temp2 = 0.0;
          for (i=1;i<=N;i++)
            Temp2 += us[i]*ut[i];
          for (i=1;i<=N;i++)
            for (j=i;j<=N;j++)
              H[i][j] += uy[i]*uy[j]/Temp1 -
                         ut[i]*ut[j]/Temp2;
        }

      }

    }


    //  -------------------------------------------------------------------

    void  BFGSMin::Choose_Lambda ( rvector X, rvector S,
                                   realtype & Lambda0 )  {
    int      i;
    realtype SS;

      for (i=1;i<=N;i++)
        if  ((S[i]!=0.0) && (!Freese[i]))  {
          SS = X[i] + Lambda0*S[i];
          if (TL)  {
            if (SS>TL[i])  Lambda0 = (TL[i]-X[i])/S[i]/(1.0+MachEps);
          }
          if (LL)  {
            if (SS<LL[i])  Lambda0 = (LL[i]-X[i])/S[i]/(1.0+MachEps);
          }
        } else if (Freese[i])  S[i] = 0.0;

    }


    //  -------------------------------------------------------------------

    void  BFGSMin::Stop()  {
      TermCode = BFGS_Stopped;
    }


    //  -------------------------------------------------------------------

    void  BFGSMin::LineSearch (  rvector   px0,   rvector   G,
                                 rvector   P,     realtype  pFunc,
                                 int  & RetCode,  bool & MaxTkn )  {

    //  A6.3.1   :   Linear  Search

    int      i;
    realtype Alpha, NewtLn, S, InitSp, RelLng, MinLam;
    realtype Lambda,LamTem, LamPre, FplPre, A, B;
    realtype Disc,  B1, B2, Lambda0;

      LamPre = 1.0;  // only to keep compiler happy
      FplPre = 1.0;  // only to keep compiler happy

      MaxTkn  = false;
      RetCode = 2;
      Alpha   = 1.0e-4;
      NewtLn  = 0.0;
      for (i=1;i<=N;i++)  {   //  calculate Newtonian step along the -gradient
        B1      = Sx[i]*P[i];
        NewtLn += B1*B1;
      }
      NewtLn = sqrt(NewtLn);

      if (NewtLn>MxStep)  {   //  restrict Newtonian step to MxStep
        S = MxStep/NewtLn;
        for (i=1;i<=N;i++)
          P[i] *= S;
        NewtLn = MxStep;
      }

      InitSp  = 0.0;
      RelLng  = 0.0;
      Lambda0 = 1.0;
      Choose_Lambda ( px0,P,Lambda0 );
      for (i=1;i<=N;i++)  {
        InitSp += G[i]*P[i];
        B1      = fabs(px0[i]);
        B2      = 1.0/Sx[i];
        if (B1>B2)  S = B1;
              else  S = B2;
        S       = fabs(P[i])/S;
        if (S>RelLng)  RelLng = S;
      }
      InitSp *= Lambda0;

      MinLam = StpEps/RelLng;
      Lambda = Lambda0;
      do {
        for (i=1;i<=N;i++)
          XPlus[i] = px0[i] + Lambda*P[i];

        MinFunc1 ( XPlus,FPlus );
        if (TermCode!=BFGS_NoTermination)  return;
        if (FPlus<=pFunc+Alpha*Lambda*InitSp)  {
          RetCode = 0;
          MaxTkn  = (Lambda==Lambda0) && (NewtLn>0.99*MxStep);
        } else if (Lambda<MinLam)  {
          RetCode = 1;
          for (i=1;i<=N;i++)
            XPlus[i] = px0[i];
        } else if (Lambda==Lambda0)  {
          LamTem = -InitSp/(2.0*(FPlus-pFunc-InitSp));
          LamTem = LamTem*Lambda;
          LamPre = Lambda;
          FplPre = FPlus;
          if (LamTem>0.1*Lambda)  Lambda = LamTem;
          else  {
            Lambda *= 0.1;
            Lambda0 = Lambda;
          }
          if (Lambda>Lambda0)  {
            Lambda  = Lambda0;
            RetCode = 0;
            for (i=1;i<=N;i++)
              XPlus[i] = px0[i] + Lambda*P[i];
          }
        } else  {
          B1 = FPlus  - pFunc - Lambda*InitSp;
          B2 = FplPre - pFunc - LamPre*InitSp;
          A  = ( B1/(Lambda*Lambda) - B2/(LamPre*LamPre) ) /
               ( Lambda - LamPre );
          B  = ( -LamPre*B1/(Lambda*Lambda) +
                 Lambda*B2/(LamPre*LamPre) ) /
               ( Lambda - LamPre );
          Disc = B*B - 3.0*A*InitSp;
          if (A==0.0)  LamTem = -InitSp/(2.0*B);
                 else  LamTem = (-B+sqrt(RMax(Disc,0.0)))/(3.0*A);
          B1 = 0.5*Lambda;
          if (B1<LamTem)  LamTem = B1;
          LamPre = Lambda;
          FplPre = FPlus;
          if (LamTem>0.1*Lambda)  Lambda = LamTem;
          else  {
            Lambda *= 0.1;
            Lambda0 = Lambda;
          }
          if (Lambda>Lambda0)  {
            Lambda  = Lambda0;
            RetCode = 0;
            for (i=1;i<=N;i++)
              XPlus[i] = px0[i] + Lambda*P[i];
          }
        }

      } while (RetCode>=2);

      TakenLambda = Lambda;

    }


    //  -------------------------------------------------------------------

    void  BFGSMin::GetMemory()  {
      if (N!=NAlloc)  {
        FreeMemory();
        GetMatrixMemory ( Hsn   , N,N, 1,1 );
        GetVectorMemory ( GPlus , N, 1 );
        GetVectorMemory ( GradX , N, 1 );
        GetVectorMemory ( HDiag , N, 1 );
        GetVectorMemory ( SN    , N, 1 );
        GetVectorMemory ( Sx    , N, 1 );
        GetVectorMemory ( XPlus , N, 1 );
        GetVectorMemory ( XOpt  , N, 1 );
        GetVectorMemory ( Freese, N, 1 );
        if (CalcHess)  {
          GetVectorMemory ( StepSize , N, 1 );
          GetVectorMemory ( FNeighbor, N, 1 );
        } else  {
          GetVectorMemory ( us       , N, 1 );
          GetVectorMemory ( uy       , N, 1 );
          GetVectorMemory ( ut       , N, 1 );
        }
        NAlloc = N;
      }
    }

    void  BFGSMin::FreeMemory()  {
      if (NAlloc>0)  {
        FreeVectorMemory ( us       , 1 );
        FreeVectorMemory ( uy       , 1 );
        FreeVectorMemory ( ut       , 1 );
        FreeVectorMemory ( Freese   , 1 );
        FreeVectorMemory ( StepSize , 1 );
        FreeVectorMemory ( FNeighbor, 1 );
        FreeVectorMemory ( XOpt     , 1 );
        FreeVectorMemory ( XPlus    , 1 );
        FreeVectorMemory ( Sx       , 1 );
        FreeVectorMemory ( SN       , 1 );
        FreeVectorMemory ( HDiag    , 1 );
        FreeVectorMemory ( GradX    , 1 );
        FreeVectorMemory ( GPlus    , 1 );
        FreeMatrixMemory ( Hsn      , NAlloc, 1,1 );
      }
      NAlloc = 0;
    }


    //  -------------------------------------------------------------------

    void  BFGSMin::Relax()  {
    int i;
      if (FPlus>FOpt)  {
        for (i=1;i<=N;i++)
          XPlus[i] = XOpt[i];
        FPlus = FOpt;
      } else  {
        for (i=1;i<=N;i++)
          XOpt[i] = XPlus[i];
        FOpt = FPlus;
      }
    }

    void  BFGSMin::CopyPlus ( rvector x0 )  {
    int i;
      for (i=1;i<=N;i++)  {
        x0   [i] = XPlus[i];
        GradX[i] = GPlus[i];
      }
      Func = FPlus;
    }


    //  -------------------------------------------------------------------

    void  BFGSMin::BFGS_Driver (  int        MinN,
                                  rvector    x0,
                                  rvector    TypX,
                                  realtype & FuncValue,
                                  int      & TerminationCode,
                                  int        Digits,
                                  int        ItnLmt,
                                  realtype   TypF,
                                  realtype   GrdTol,
                                  realtype   StpTol,
                                  realtype   MaxStp,
                                  bool       Hess,
                                  rvector    LowLimit,
                                  rvector    TopLimit )  {

    //  D6.1.1   :  Unconstrained Minimization Driver

    int   i,RetCode;
    int   ItnCnt;
    bool  MaxTkn;

      TL       = TopLimit;
      LL       = LowLimit;
      ForDiff  = true;
      N        = MinN;
      CalcHess = Hess;

      ModF     = false;

      GetMemory();

      UMInCk ( x0,TypX,Digits,TypF,
               GrdTol,StpTol,MaxStp,
               ItnLmt );
      if (TermCode!=BFGS_NoTermination)  {
        FreeMemory();
        FuncValue       = Func;
        TerminationCode = TermCode;
        return;
      }

      ItnCnt = 0;

      MinFunc1 ( x0,Func );
      if (TermCode!=BFGS_NoTermination)  {
        FreeMemory();
        FuncValue       = Func;
        TerminationCode = TermCode;
        return;
      }
      FOpt  = Func;
      FPlus = Func;
      for (i=1;i<=N;i++)  {
        XOpt [i] = x0[i];
        XPlus[i] = x0[i];
      }
      ModF = true;
      Gradient ( x0,GradX,Func );
      Print    ( ItnCnt,x0,GradX,Func );
      for (i=1;i<=N;i++)
        GPlus[i] = GradX[i];
      if (TermCode!=BFGS_NoTermination)  {
        Relax     ();
        CopyPlus  ( x0 );
        FreeMemory();
        FuncValue       = Func;
        TerminationCode = TermCode;
        return;
      }

      UMStop0 ( x0,GradX );
      if  (TermCode!=BFGS_NoTermination)  {
        FreeMemory();
        FuncValue       = Func;
        TerminationCode = TermCode;
        return;
      }

      if (!CalcHess)  InitHessUnFac ( Func,Hsn );

      RetCode = 0;
      while (TermCode==BFGS_NoTermination)  {
        ItnCnt++;
        if (RetCode>=0)  {
          if (CalcHess)  {
            FDHessF ( Func,x0 );
            if (TermCode!=BFGS_NoTermination)  {
              Relax     ();
              CopyPlus  ( x0 );
              FreeMemory();
              FuncValue       = Func;
              TerminationCode = TermCode;
              return;
            }
          }
          MdHess ( Hsn,HDiag );
        }
        ChSolve    ( N,Hsn,GradX,SN );
        LineSearch ( x0,GradX,SN,Func,RetCode,MaxTkn );
        if ((RetCode==1) && ForDiff)  {
          RetCode = -1;
          ForDiff = false;
        } else
          Relax();
        if (TermCode!=BFGS_NoTermination)  {
          Relax     ();
          CopyPlus  ( x0 );
          FreeMemory();
          FuncValue       = Func;
          TerminationCode = TermCode;
          return;
        } else
          Gradient ( XPlus,GPlus,FPlus );
        if (TermCode!=BFGS_NoTermination)  {
          Relax     ();
          CopyPlus  ( x0 );
          FreeMemory();
          FuncValue       = Func;
          TerminationCode = TermCode;
          return;
        }
        if (RetCode>=0)  {
          UMStop ( x0,GPlus,RetCode,ItnCnt,MaxTkn );
          if ((!CalcHess) && (TermCode==BFGS_NoTermination))
            BFGSUnFac ( x0,XPlus,GradX,GPlus,false,HDiag,Hsn );
        }
        CopyPlus ( x0 );
        Print    ( ItnCnt, x0,GradX,Func );
      }

      Relax     ();
      FreeMemory();
      FuncValue       = Func;
      TerminationCode = TermCode;

    }

  }  // namespace math

}  // namespace mmdb

