/* Copyright (C) 2005-2017 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>

#include <libscuff.h>
#include <libIncField.h>
#include <libSpherical.h>
#include "GLMTSolver.h"
using namespace scuff;

#define II cdouble(0.0,1.0)

namespace scuff{
HMatrix *GetSCRMatrix(char *BSMesh, double R, int NumPoints,
                      GTransformation *GT1=0, GTransformation *GT2=0);
void GetNMatrices(double nHat[3], double X[3], double XTorque[3],
                  double NMatrix[NUMPFT][3][3], 
                  bool *NeedQuantity=0);
double HVMVP(cdouble V1[3], double M[3][3], cdouble V2[3]);
               }

/***************************************************************/
/* identify which region of a spherically-symmetric geometry   */
/* we're in                                                    */
/* region 0 is outermost (exterior), region NumRegions-1       */
/* is innermost                                                */
/***************************************************************/
int GetRegion(GLMTSolver *GLMTS, double r)
{
  int NumRegions = GLMTS->NumRegions;
  double *Radii  = GLMTS->Radii;

  for(int nr=1; nr<NumRegions; nr++)
   if (r>Radii[nr])
    return nr-1;

  return NumRegions-1;
}

int GetRegion(GLMTSolver *GLMTS, double X[3])
{ return GetRegion(GLMTS, VecNorm(X)); }

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble GetEps(GLMTSolver *GLMTS, double Omega, double X[3])
{
  return GLMTS->MPs[ GetRegion(GLMTS, X) ] -> GetEps(Omega);
}

/***************************************************************/
/* layout of coefficients:   ***********************************/
/*  AHat_{0,Alpha=1}          **********************************/
/*  A   _{1,Alpha=1}          **********************************/
/*  AHat_{1,Alpha=1}          **********************************/
/*  A   _{2,Alpha=1}          **********************************/
/*  AHat_{2,Alpha=1}          **********************************/
/*  ...                       **********************************/
/*  A   _{NumRegions-1, Alpha=1}                               */
/*  AHat_{0,Alpha=2}                                           */
/*  ...                                                        */
/*  here Alpha is the spherical-wave index:                    */
/*   Alpha | (l,m)                                             */
/*       0 | (0,0)                                             */
/*       1 | (1,-1)                                            */
/*       2 | (1, 0)                                            */
/*       3 | (1,+1)                                            */
/*       4 | (2,-2)                                            */
/*     ... | ...                                               */
/***************************************************************/
#define COEFF_AREG 0
#define COEFF_AHAT 1
#define COEFF_BREG 2
#define COEFF_BHAT 3
int GetCoeffIndex(GLMTSolver *GLMTS, int Alpha,
                  int WhichRegion, int RegHat)
{
  int NumRegions     = GLMTS->NumRegions;
  int CoeffsPerAlpha = 2*(NumRegions-1);
  int LMax           = GLMTS->LMax;
  int NAlpha         = (LMax+1)*(LMax+1);
  
  int Offset = RegHat ? 1 : 0;
  int Coeff = -1 + 2*WhichRegion + Offset;
  if (Coeff<0 || Coeff>CoeffsPerAlpha)
   return -1;

  return Alpha*CoeffsPerAlpha + Coeff;
}

cdouble GetCoefficient(GLMTSolver *GLMTS, int Alpha, int WhichRegion,
                       int WhichCoefficient)
{
  HVector *CVector;
  int RegHat;
  switch(WhichCoefficient)
   { case COEFF_AREG: CVector = GLMTS->AVector; RegHat=0; break;
     case COEFF_AHAT: CVector = GLMTS->AVector; RegHat=1; break;
     case COEFF_BREG: CVector = GLMTS->BVector; RegHat=0; break;
     case COEFF_BHAT: CVector = GLMTS->BVector; RegHat=1; break;
   };

  int Index = GetCoeffIndex(GLMTS, Alpha, WhichRegion, RegHat);
  if (Index<0 || Index>=CVector->N) return 0.0;
  return CVector->ZV[Index];
  
}

/***************************************************************/
/* fill in the PVector and QVector fields in GLMTS with the    */
/* M- and N-type spherical-wave expansion coefficients for     */
/* the incident field described by IF                          */
/* (only works for IF = spherical wave or plane wave traveling */
/*  in positive z-direction)                                   */
/***************************************************************/
#define IPOW(n) ( (n)==0 ? 1.0 : (n)==1 ? II : (n)==2 ? -1.0 : -1.0*II)
void GetRHSExpansion(GLMTSolver *GLMTS, IncField *IF)
{
  int LMax         = GLMTS->LMax;

  HVector *PVector=GLMTS->PVector;
  HVector *QVector=GLMTS->QVector;

  PVector->Zero();
  QVector->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SphericalWave *SW = dynamic_cast<SphericalWave *>(IF);
  if (SW)
   { int L = SW->L;
     int M = SW->M;
     int Alpha = LM2ALPHA(L,M);
     if (SW->Type == SW_MAGNETIC)
      PVector->SetEntry(Alpha,1.0);
     else
      QVector->SetEntry(Alpha,1.0);
     return;
   };
 
  /***************************************************************/
  /* use known expansion for z-traveling plane wave              */
  /***************************************************************/
  PlaneWave *PW = dynamic_cast<PlaneWave *>(IF);
  double zHat[3]={0.0, 0.0, 1.0};
  if (PW && VecEqualFloat(PW->nHat, zHat))
   {
     cdouble Ex = PW->E0[0];
     cdouble Ey = PW->E0[1];
     double EMag = sqrt( abs(Ex)*abs(Ex) + abs(Ey)*abs(Ey) );
     bool RCP = EqualFloat( real(Ex),      imag(Ey) );
     bool LCP = EqualFloat( real(Ex), -1.0*imag(Ey) );
     bool LP  = ( imag(Ey)==0.0 );
     double Factor = (LP ? 0.5 : 1.0/M_SQRT2)*EMag;
     for(int L=1; L<=LMax; L++)
      { 
        cdouble Pl   = Factor*IPOW( (L%4) )*sqrt(4.0*M_PI*(2.0*L+1.0));
        if (RCP || LP)
         { 
           int Alpha = LM2ALPHA(L,+1);
           PVector->AddEntry(Alpha, Pl);
           QVector->AddEntry(Alpha, -1.0*II*Pl);
         };
        if (LCP || LP)
         { 
           int Alpha = LM2ALPHA(L,-1);
           PVector->AddEntry(Alpha, Pl);
           QVector->AddEntry(Alpha, +1.0*II*Pl);
         };
      };
   };
}

/***************************************************************/
/* get the total, scattered, and incident E, H fields at X *****/
/***************************************************************/
void GetFields(GLMTSolver *GLMTS, double Omega, IncField *IF,
               double X[3], 
               cdouble EHTot[6], cdouble EHScat[6], cdouble EHInc[6])
{
  int NumRegions     = GLMTS->NumRegions;
  int LMax           = GLMTS->LMax;
  cdouble *MRegArray = GLMTS->MRegArray;
  cdouble *MOutArray = GLMTS->MOutArray;
  cdouble *NRegArray = GLMTS->NRegArray;
  cdouble *NOutArray = GLMTS->NOutArray;
  HVector *PVector   = GLMTS->PVector;
  HVector *QVector   = GLMTS->QVector;
   double *Workspace = GLMTS->Workspace;

  /***************************************************************/
  /* figure out what region we're in *****************************/
  /***************************************************************/
  double r, Theta, Phi;
  CoordinateC2S(X, &r, &Theta, &Phi);

  int nr = GetRegion(GLMTS, r);
  cdouble ZRel, k = Omega*GLMTS->MPs[nr]->GetRefractiveIndex(Omega, &ZRel);
  
  bool Innermost = (nr==(NumRegions-1));
  bool Outermost = (nr==0);

  /***************************************************************/
  /* fetch M, N functions for this region ************************/
  /***************************************************************/
  GetMNlmArray(LMax, k, r, Theta, Phi, LS_REGULAR,
               MRegArray, NRegArray, Workspace);

  if (!Innermost)
   GetMNlmArray(LMax, k, r, Theta, Phi, LS_OUTGOING,
                MOutArray, NOutArray, Workspace);

  /***************************************************************/
  /* sum contributions of all waves to fields in this region     */
  /***************************************************************/
  cdouble EHScatS[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  cdouble EHIncS[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for(int L=1, Alpha=1; L<=LMax; L++)
   for(int M=-L; M<=L; M++, Alpha++)
    {  
      cdouble AReg = GetCoefficient(GLMTS, Alpha, nr, COEFF_AREG);
      cdouble AOut = GetCoefficient(GLMTS, Alpha, nr, COEFF_AHAT);
      cdouble BReg = GetCoefficient(GLMTS, Alpha, nr, COEFF_BREG);
      cdouble BOut = GetCoefficient(GLMTS, Alpha, nr, COEFF_BHAT);

      if (Outermost)
       { cdouble P = PVector->GetEntry(Alpha);
         cdouble Q = QVector->GetEntry(Alpha);
         VecPlusEquals(EHIncS+0,         P, MRegArray + 3*Alpha);
         VecPlusEquals(EHIncS+0,         Q, NRegArray + 3*Alpha);
         VecPlusEquals(EHIncS+3,         Q, MRegArray + 3*Alpha);
         VecPlusEquals(EHIncS+3,    -1.0*P, NRegArray + 3*Alpha);
       }
      else
       { VecPlusEquals(EHScatS+0,      AReg, MRegArray + 3*Alpha);
         VecPlusEquals(EHScatS+0,      BReg, NRegArray + 3*Alpha);
         VecPlusEquals(EHScatS+3,      BReg, MRegArray + 3*Alpha);
         VecPlusEquals(EHScatS+3, -1.0*AReg, NRegArray + 3*Alpha);
       }
 
      if (!Innermost)
       { VecPlusEquals(EHScatS+0,      AOut, MOutArray + 3*Alpha);
         VecPlusEquals(EHScatS+0,      BOut, NOutArray + 3*Alpha);
         VecPlusEquals(EHScatS+3,      BOut, MOutArray + 3*Alpha);
         VecPlusEquals(EHScatS+3, -1.0*AOut, NOutArray + 3*Alpha);
       };
    };

  VectorS2C(Theta, Phi, EHScatS+0, EHScat+0);
  VectorS2C(Theta, Phi, EHScatS+3, EHScat+3);
  VecScale(EHScat+3, 1.0/(ZVAC*ZRel), 3);

  VectorS2C(Theta, Phi, EHIncS+0, EHInc+0);
  VectorS2C(Theta, Phi, EHIncS+3, EHInc+3);
  VecScale(EHInc+3, 1.0/ZVAC, 3);

  memcpy(EHTot, EHScat, 6*sizeof(cdouble));
  if (Outermost)
   VecPlusEquals(EHTot, 1.0, EHInc, 6);

}

void GetFields(GLMTSolver *GLMTS, double Omega, IncField *IF,
               double X[3], cdouble EHTot[6])
{
  cdouble EHScat[6], EHInc[6];
  GetFields(GLMTS, Omega, IF, X, EHTot, EHScat, EHInc);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void TestIncFields(GLMTSolver *GLMTS, double Omega, IncField *IF,
                   double X[3], cdouble EHExact[6], cdouble EHExpand[6])
{
  cdouble *MRegArray = GLMTS->MRegArray;
  cdouble *NRegArray = GLMTS->MOutArray;
   double *Workspace = GLMTS->Workspace;

  int LMax=GLMTS->LMax;
  int NAlpha = (LMax+1)*(LMax+1);

  /***************************************************************/
  /* figure out what region we're in *****************************/
  /***************************************************************/
  double r, Theta, Phi;
  CoordinateC2S(X, &r, &Theta, &Phi);

  GetMNlmArray(LMax, Omega, r, Theta, Phi, LS_REGULAR,
               MRegArray, NRegArray, Workspace);

  HVector *PVector = GLMTS->PVector;
  HVector *QVector = GLMTS->QVector;
  GetRHSExpansion(GLMTS, IF);

  cdouble EHS[6];
  memset(EHS,0,6*sizeof(cdouble));
  for(int L=1, Alpha=1; L<=LMax; L++)
   for(int M=-L; M<=L; M++, Alpha++)
    {  
      cdouble P = PVector->GetEntry(Alpha);
      cdouble Q = QVector->GetEntry(Alpha);
      VecPlusEquals(EHS+0, P, MRegArray + 3*Alpha);
      VecPlusEquals(EHS+0, Q, NRegArray + 3*Alpha);
      VecPlusEquals(EHS+3, Q, MRegArray + 3*Alpha);
      VecPlusEquals(EHS+3, -1.0*P, NRegArray + 3*Alpha);
    };

  VectorS2C(Theta, Phi, EHS+0, EHExpand+0);
  VectorS2C(Theta, Phi, EHS+3, EHExpand+3);
  VecScale(EHExpand+3, 1.0/ZVAC, 3);

  IF->GetFields(X, EHExact);

}

/***************************************************************/
/* displaced surface-integral power, force, torque             */
/***************************************************************/
void GetDSIPFT(GLMTSolver *GLMTS, double Omega, IncField *IF,
               double PFT[NUMPFT], double DSIRadius, int DSIPoints)
{
  Log("Computing DSIPFT: (R,NPts)=(%e,%i)",DSIRadius,DSIPoints);

  /***************************************************************/
  /* get cubature-rule matrix ************************************/
  /***************************************************************/
  HMatrix *SCRMatrix = GetSCRMatrix(0, DSIRadius, DSIPoints, 0, 0);

  /***************************************************************/
  /* we assume that all cubature points lie in the same region   */
  /* of the scuff geometry, so we use the first point in the rule*/
  /* to determine which region that is and look up its eps/mu    */
  /***************************************************************/
  double EpsAbs = TENTHIRDS / ZVAC;
  double  MuAbs = TENTHIRDS * ZVAC;

  double XTorque[3] = {0.0, 0.0, 0.0};

  /***************************************************************/
  /* loop over points in the cubature rule                       */
  /***************************************************************/
  Log(" Doing DSI calculation...");
  memset(PFT, 0, NUMPFT*sizeof(double));
  for(int nr=0; nr<SCRMatrix->NR; nr++)
   { 
     LogPercent(nr, SCRMatrix->NR, 10);

     double w, X[3], nHat[3];
     SCRMatrix->GetEntriesD(nr, "0:2", X);
     SCRMatrix->GetEntriesD(nr, "3:5", nHat);
     w = SCRMatrix->GetEntryD(nr, 6);

     double NMatrix[NUMPFT][3][3];
     GetNMatrices(nHat, X, XTorque, NMatrix);
     
     cdouble EHT[6], EHS[6], EHI[6];
     cdouble *ET = EHT+0, *HT=EHT+3;
     cdouble *ES = EHS+0, *HS=EHS+3;
     GetFields(GLMTS, Omega, 0, X, EHT, EHS, EHI);

     // scattered power
     double dP = +0.25 * w * (  HVMVP(ES, NMatrix[PFT_PABS], HS)
                               -HVMVP(HS, NMatrix[PFT_PABS], ES)
                             );
     PFT[PFT_PSCAT] += dP;

     // absorbed power
     dP = -0.25 * w * (  HVMVP(ET, NMatrix[PFT_PABS], HT)
                        -HVMVP(HT, NMatrix[PFT_PABS], ET)
                      );
     PFT[PFT_PABS] += dP;

     // force and torque
     double dFT[NUMPFT];
     for(int nq=2; nq<NUMPFT; nq++)
      { dFT[nq] = 0.25 * w * ( EpsAbs*HVMVP(ET, NMatrix[nq], ET)
                               +MuAbs*HVMVP(HT, NMatrix[nq], HT)
                             );
        PFT[nq] += dFT[nq];
      };

   };
  delete SCRMatrix;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SetEntry(HMatrix *M, int nr, int nc, cdouble Entry)
{ if ( nr>=0 && nr<M->NR && nc>=0 && nc<M->NC )
   M->SetEntry(nr,nc,Entry);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Solve(GLMTSolver *GLMTS, double Omega, IncField *IF)
{
  int NumRegions     = GLMTS->NumRegions;
  double *Radii      = GLMTS->Radii;
  MatProp **MPs      = GLMTS->MPs;
  HVector *AVector   = GLMTS->AVector;
  HVector *BVector   = GLMTS->BVector;
  AVector->Zero();
  BVector->Zero();

  int LMax           = GLMTS->LMax;
  int NAlpha         = (LMax+1)*(LMax+1);

  /***************************************************************/
  /* RRegMatrix[L, 1]   = RReg_L(k_0 * r1)                       */
  /* RRegMatrix[L, 2]   = RReg_L(k_1 * r1)                       */
  /* RRegMatrix[L, 3]   = RReg_L(k_1 * r2)                       */
  /* RRegMatrix[L, 4]   = RReg_L(k_2 * r2)                       */
  /* RRegMatrix[L, 5]   = RReg_L(k_2 * r3)                       */
  /* ...                                                         */
  /* RRegMatrix[L, 2*n] = RReg_L(k_n * rn)                       */
  /***************************************************************/
  HMatrix *RRegMatrix    = new HMatrix(LMax+2, 2*NumRegions, LHM_COMPLEX);
  HMatrix *RBarRegMatrix = new HMatrix(LMax+2, 2*NumRegions, LHM_COMPLEX);
  HMatrix *ROutMatrix    = new HMatrix(LMax+2, 2*NumRegions, LHM_COMPLEX);
  HMatrix *RBarOutMatrix = new HMatrix(LMax+2, 2*NumRegions, LHM_COMPLEX);
  HVector *ZRegion       = new HVector(NumRegions, LHM_COMPLEX);
  double *Workspace      = GLMTS->Workspace;
  for(int nr=0; nr<NumRegions; nr++)
   { 
     cdouble ZRel, k = Omega * MPs[nr]->GetRefractiveIndex(Omega, &ZRel);

     ZRegion->SetEntry(nr, ZRel);

     for(int p=0; p<2; p++)
      { 
        if ( nr==0 && p==0 ) continue;
        if ( nr==(NumRegions-1) && p==1 ) continue;
        double r = Radii[nr+p];

        cdouble *RReg    = (cdouble *)RRegMatrix->GetColumnPointer(2*nr+p);
        cdouble *RBarReg = (cdouble *)RBarRegMatrix->GetColumnPointer(2*nr+p);
        cdouble *ROut    = (cdouble *)ROutMatrix->GetColumnPointer(2*nr+p);
        cdouble *RBarOut = (cdouble *)RBarOutMatrix->GetColumnPointer(2*nr+p);

        GetRadialFunctions(LMax, k, r, LS_REGULAR,  RReg, RBarReg, Workspace);
        GetRadialFunctions(LMax, k, r, LS_OUTGOING, ROut, RBarOut, Workspace);
        for(int L=1; L<=LMax; L++)
         { RBarReg[L] = RBarReg[L]/k + RReg[L]/(k*r);
           RBarOut[L] = RBarOut[L]/k + ROut[L]/(k*r);
         };
      };
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HVector *PVector = GLMTS->PVector;
  HVector *QVector = GLMTS->QVector;
  GetRHSExpansion(GLMTS, IF);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int Dim = 2*(NumRegions-1);
  HMatrix *AMatrix    = new HMatrix(Dim, Dim, LHM_COMPLEX);
  HMatrix *BMatrix    = new HMatrix(Dim, Dim, LHM_COMPLEX);
  HVector *PartialRHS = new HVector(Dim, LHM_COMPLEX);
  for(int L=1; L<=LMax; L++)
   { 
     /*--------------------------------------------------------------*/
     /*- form and factorize matrices for this L ---------------------*/
     /*--------------------------------------------------------------*/
     AMatrix->Zero();
     BMatrix->Zero();
     for(int nr=0; nr<NumRegions; nr++)
      for(int p=0; p<2; p++)
       { 
         if (nr==0 && p==0) continue;
         if (nr==(NumRegions-1) && p==1) continue;

         cdouble ZRel = ZRegion->GetEntry(nr);
         int ncReg = 2*nr - 1,     ncHat = 2*nr;

         int nrE = 2*(nr-1+p) + 0, nrH   = 2*(nr-1+p) + 1;
         int nR = 2*nr + p;
         double Sign = (p==1) ? -1.0 : 1.0;

         // E-field match at region nr+p, nr+p+1 boundary
         SetEntry(AMatrix, nrE, ncReg, Sign*   RRegMatrix->GetEntry(L,nR) );
         SetEntry(BMatrix, nrE, ncReg, Sign*RBarRegMatrix->GetEntry(L,nR) );
         SetEntry(AMatrix, nrE, ncHat, Sign*   ROutMatrix->GetEntry(L,nR) );
         SetEntry(BMatrix, nrE, ncHat, Sign*RBarOutMatrix->GetEntry(L,nR) );

         // H-field match at region nr+p, nr+p boundary
         SetEntry(AMatrix, nrH, ncReg, Sign*RBarRegMatrix->GetEntry(L,nR)/ZRel );
         SetEntry(BMatrix, nrH, ncReg, Sign*   RRegMatrix->GetEntry(L,nR)/ZRel );
         SetEntry(AMatrix, nrH, ncHat, Sign*RBarOutMatrix->GetEntry(L,nR)/ZRel );
         SetEntry(BMatrix, nrH, ncHat, Sign*   ROutMatrix->GetEntry(L,nR)/ZRel );
      }; // for(int nr=0; nr<NumRegions; nr++)
     AMatrix->LUFactorize();
     BMatrix->LUFactorize();

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     cdouble ZOut=ZRegion->GetEntry(0);
     for(int M=-L; M<=L; M++)
      { 
        int Alpha = LM2ALPHA(L, M);

        cdouble P=PVector->GetEntry(Alpha);
        if (abs(P)>0.0)
         { PartialRHS->Zero();
           PartialRHS->SetEntry(0, P*RRegMatrix->GetEntry(L,1));
           PartialRHS->SetEntry(1, P*RBarRegMatrix->GetEntry(L,1)/ZOut);
           AMatrix->LUSolve(PartialRHS);
           for(int d=0; d<Dim; d++)
            AVector->SetEntry(Alpha*Dim + d, PartialRHS->GetEntry(d));
         };

        cdouble Q=QVector->GetEntry(Alpha);
        if (abs(Q)>0.0)
         { PartialRHS->Zero();
           PartialRHS->SetEntry(0, Q*RBarRegMatrix->GetEntry(L,1));
           PartialRHS->SetEntry(1, Q*RRegMatrix->GetEntry(L,1)/ZOut);
           BMatrix->LUSolve(PartialRHS);
           for(int d=0; d<Dim; d++)
            BVector->SetEntry(Alpha*Dim + d, PartialRHS->GetEntry(d));
         };

      }; // for(int M=-L; M<=L; M++)

   }; // for(int L=1, Alpha=1; L<=LMax; L++)

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  delete AMatrix;
  delete BMatrix;
  delete PartialRHS;
  delete RRegMatrix;
  delete RBarRegMatrix;
  delete ROutMatrix;
  delete RBarOutMatrix;
  delete ZRegion;
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
GLMTSolver *CreateGLMTSolver(char *GLMTFile, int LMax)
{
  FILE *f=fopen(GLMTFile,"r");
  if (!f) ErrExit("could not open %s",GLMTFile);

  #define MAXREGIONS 10
  double Radii[MAXREGIONS];
  MatProp *MPs[MAXREGIONS];
 
  Radii[0]       = 1.0e10;
  MPs[0]         = new MatProp(MP_VACUUM);
  int NumRegions = 1;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  char Line[100];
  int LineNum=0;
  while(fgets(Line, 100, f))
   { 
     LineNum++;

     if (strlen(Line)==0) continue;

     Line[strlen(Line)-1]=0;
     char *Tokens[2];
     int NumTokens = Tokenize(Line, Tokens, 2);
     if (NumTokens==0 || Tokens[0][0]=='#') 
      continue;

     double R;
     char *MP;
     if ( NumTokens!=2 )
      ErrExit("%s:%i: syntax error",GLMTFile,LineNum);

     if (!strcasecmp(Tokens[0],"MEDIUM"))
      { MPs[0] = new MatProp(Tokens[1]);
        if (MPs[0]->ErrMsg) 
         ErrExit("%s:%i: %s: %s",GLMTFile,LineNum,Tokens[1],MPs[NumRegions]->ErrMsg);
        printf("Exterior region: material=%s\n",MPs[0]->Name);
      }
     else
      {
        if (NumRegions==MAXREGIONS)
         ErrExit("%s:%i: too many regions",GLMTFile,LineNum);
        if ( 1!=sscanf(Tokens[0],"%le",Radii + NumRegions) )
         ErrExit("%s:%i: invalid radius %s",GLMTFile,LineNum,Tokens[0]);
        MPs[NumRegions] = new MatProp(Tokens[1]);
        if (MPs[NumRegions]->ErrMsg) 
         ErrExit("%s:%i: %s: %s",GLMTFile,LineNum,Tokens[1],MPs[NumRegions]->ErrMsg);
        printf("Region %i: R=%e, material=%s\n",
               NumRegions,Radii[NumRegions], MPs[NumRegions]->Name);
        NumRegions++;
      };
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  GLMTSolver *GLMTS   = (GLMTSolver *)mallocEC(sizeof(GLMTSolver));
  GLMTS->LMax       = LMax;

  GLMTS->NumRegions = NumRegions;
  GLMTS->Radii      = (double *)memdup(Radii, NumRegions*sizeof(double));
  GLMTS->MPs        = (MatProp **)memdup(MPs, NumRegions*sizeof(MatProp *));
  
  int NAlpha        = (LMax+1)*(LMax+1);
  int NCoeff        = 2*(NumRegions-1)*NAlpha;
  GLMTS->AVector    = new HVector(NCoeff, LHM_COMPLEX);
  GLMTS->BVector    = new HVector(NCoeff, LHM_COMPLEX);
  GLMTS->PVector    = new HVector(NCoeff, LHM_COMPLEX);
  GLMTS->QVector    = new HVector(NCoeff, LHM_COMPLEX);
   
  GLMTS->MRegArray  = (cdouble *)mallocEC(3*NAlpha*sizeof(cdouble));
  GLMTS->MOutArray  = (cdouble *)mallocEC(3*NAlpha*sizeof(cdouble));
  GLMTS->NRegArray  = (cdouble *)mallocEC(3*NAlpha*sizeof(cdouble));
  GLMTS->NOutArray  = (cdouble *)mallocEC(3*NAlpha*sizeof(cdouble));
  GLMTS->Workspace  = (double *)mallocEC( 4*(LMax+2)*sizeof(double));

  return GLMTS;
  

