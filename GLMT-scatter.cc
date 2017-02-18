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

#include <libIncField.h>
#include "GLMTSolver.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXPW    10    // max number of plane waves
#define MAXSW    10    // max number of spherical waves
#define MAXEPF   10    // max number of evaluation-point files

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteFields(GLMTSolver *GLMTS, double Omega, IncField *IF,
                 char *EPFile, char *FileBase)
{
  HMatrix *XMatrix=new HMatrix(EPFile);
  FILE *f=vfopen("%s.%s.GLMTFields","w",FileBase,GetFileBase(EPFile));
  fprintf(f,"#1-3   x,y,z\n");
  fprintf(f,"#4     Omega\n");
  fprintf(f,"#5,6   real, imag Eps\n");
  fprintf(f,"#7,8   real, imag Ex total\n");
  fprintf(f,"#9,10  real, imag Ey total \n");
  fprintf(f,"#11,12 real, imag Ez total \n");
  fprintf(f,"#13,14 real, imag Hx total \n");
  fprintf(f,"#15,16 real, imag Hy total \n");
  fprintf(f,"#17,18 real, imag Hz total \n");
  fprintf(f,"#19,20 real, imag Jx total \n");
  fprintf(f,"#21,22 real, imag Jy total \n");
  fprintf(f,"#23,24 real, imag Jz total \n");  
  fprintf(f,"#25,26 real, imag Ex incident\n");
  fprintf(f,"#27,28 real, imag Ey incident\n");
  fprintf(f,"#29,30 real, imag Ez incident\n");
  fprintf(f,"#31,32 real, imag Hx incident\n");
  fprintf(f,"#33,34 real, imag Hy incident\n");
  fprintf(f,"#35,36 real, imag Hz incident\n");

  for(int nx=0; nx<XMatrix->NR; nx++)
   { 
     double X[3];
     XMatrix->GetEntriesD(nx,":",X);
 
     cdouble Eps=GetEps(GLMTS, Omega, X);
     cdouble Factor = -1.0*II*ZVAC*Omega;

     cdouble EH[6], EHScat[6], EHInc[6], J[3];
     GetFields(GLMTS, Omega, IF, X, EH, EHScat, EHInc);

     fprintVec(f,X);
     fprintf(f,"%e %e %e ",Omega,real(Eps),imag(Eps));
     fprintVec(f,EH,6);
     EH[0] *= Factor;
     EH[1] *= Factor;
     EH[2] *= Factor;
     fprintVec(f,EH,3);
     fprintVecCR(f,EHInc,6);
   };
  fclose(f);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();
  InitializeLog(argv[0]);
  char *s;
  if ( (s=getenv("SCUFF_ABORT_ON_FPE")) && (s[0]=='1') )
   feenableexcept(FE_INVALID | FE_OVERFLOW);

  /***************************************************************/
  /* process remaining options ***********************************/
  /***************************************************************/
  char *GLMTFile=0;
  int LMax=3;
//
  double Omega=0.1;
  char *OmegaFile=0;
//
  double pwDir[3*MAXPW];             int npwDir;
  cdouble pwPol[3*MAXPW];            int npwPol;
//
  int swLMP[3*MAXSW];                int nswLMP;
//
  char *EPFiles[MAXEPF];             int nEPFiles;
  char *FileBase=0;
//
  char *PFTFile=0;
  double DSIRadius=3.0;
  int DSIPoints=302;
//
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"GLMTFile",       PA_STRING,  1, 1,       (void *)&GLMTFile,   0,       ""},
     {"LMax",           PA_INT,     1, 1,       (void *)&LMax,       0,       ""},
     {"Omega",          PA_DOUBLE,  1, 1,       (void *)&Omega,      0,       ""},
     {"OmegaFile",      PA_STRING,  1, 1,       (void *)&OmegaFile,  0,       ""},
/**/
     {"pwDirection",    PA_DOUBLE,  3, MAXPW,   (void *)pwDir,       &npwDir, "plane wave direction"},
     {"pwPolarization", PA_CDOUBLE, 3, MAXPW,   (void *)pwPol,       &npwPol, "plane wave polarization"},
/**/
     {"swLMP",          PA_INT,     3, MAXSW,   (void *)swLMP,       &nswLMP, "spherical-wave L, M, P indices"},
/**/
     {"EPFile",         PA_STRING,  1, MAXEPF,  (void *)EPFiles,     &nEPFiles, "list of evaluation points"},
/**/
     {"FileBase",       PA_STRING,  1, 1,       (void *)&FileBase, 0,         "base name of output files"},
     {"PFTFile",        PA_STRING,  1, 1,       (void *)&PFTFile, 0,          ""},
     {"DSIPoints",      PA_INT,     1, 1,       (void *)&DSIPoints,  0,       ""},
     {"DSIRadius",      PA_DOUBLE,  1, 1,       (void *)&DSIRadius,  0,       ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  GLMTSolver *GLMTS = CreateGLMTSolver(GLMTFile, LMax);
  if (!FileBase) 
   FileBase=strdup(GetFileBase(GLMTFile));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HVector *OmegaList;
  if (OmegaFile)
   OmegaList=new HVector(OmegaFile);
  else
   { OmegaList=new HVector(1);
     OmegaList->SetEntry(0,Omega);
   };

  /*******************************************************************/
  /* process incident-field-related options to construct the data    */
  /* used to generate the incident field in our scattering problem.  */
  /*******************************************************************/
  if ( npwPol != npwDir )
   ErrExit("numbers of --pwPolarization and --pwDirection options must agree");

  IncField *IF=0;
  for(int npw=0; npw<npwPol; npw++)
   { PlaneWave *PW=new PlaneWave(pwPol + 3*npw, pwDir + 3*npw);
     PW->Next=IF;
     IF=PW;
   };

  for(int nsw=0; nsw<nswLMP; nsw++)
   { int L = swLMP[3*nsw + 0];
     int M = swLMP[3*nsw + 1];
     int P = swLMP[3*nsw + 2];
     int Type = (P==0) ? SW_MAGNETIC : SW_ELECTRIC;
     SphericalWave *SW=new SphericalWave(L, M, Type);
     SW->Next=IF;
     IF=SW;
   };

  if (IF==0)
   ErrExit("you must specify an incident field");

  /***************************************************************/ 
  /***************************************************************/
  /***************************************************************/
  if (PFTFile)
   { FILE *f=fopen(PFTFile,"w");
     fprintf(f,"#1 omega\n");
     fprintf(f,"#2 PAbs\n");
     fprintf(f,"#3 PScat\n");
     fprintf(f,"#4-6 Fx, Fy, Fz\n");
     fprintf(f,"#7-9 Tx, Ty, Tz\n");
     fclose(f);
   };

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  for(int nOmega=0; nOmega<OmegaList->N; nOmega++)
   { Omega = OmegaList->GetEntryD(nOmega);
     IF->SetFrequency(Omega);

     Solve(GLMTS, Omega, IF);

     for(int nepf=0; nepf<nEPFiles; nepf++)
      WriteFields(GLMTS, Omega, IF, EPFiles[nepf], FileBase);

     if (PFTFile)
      { double PFT[NUMPFT];
        GetDSIPFT(GLMTS, Omega, IF, PFT, DSIRadius, DSIPoints);
        FILE *f=fopen(PFTFile,"a");
        fprintf(f,"%e ",Omega);
        fprintVecCR(f,PFT,NUMPFT);
        fclose(f);
      };
   };

}
