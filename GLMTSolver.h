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

#ifndef GLMT_H 
#define GLMT_H

#include <libhrutil.h>
#include <libhmat.h>
#include <libMatProp.h>
#include <libIncField.h>
#include <PFTOptions.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct GLMTSolver
{
  // info on the geometry 
  int LMax;
  int NumRegions;
  double *Radii;
  MatProp **MPs;

  // vectors of spherical-wave expansion coefficients for the
  // fields in the various regions 
  HVector *AVector, *BVector, *PVector, *QVector;
    
  // internal storage buffers
  cdouble *MRegArray, *MOutArray, *NRegArray, *NOutArray;
  double *Workspace;

} GLMTSolver;

/***************************************************************/
/* create a GLMTSolver structure from a .GLMT file *************/
/***************************************************************/
GLMTSolver *CreateGLMTSolver(char *GLMTFile, int LMax);

/***************************************************************/
/* solve the scattering problem for a given frequency and      */
/* incident field                                              */
/***************************************************************/
void Solve(GLMTSolver *GLMTS, double Omega, IncField *IF);

/***************************************************************/
/* get components of E, H fields at arbitrary points ***********/
/***************************************************************/
void GetFields(GLMTSolver *GLMTS, double Omega, IncField *IF,
               double X[3], cdouble EHTot[6], 
               cdouble EHScat[6], cdouble EHInc[6]);

void GetFields(GLMTSolver *GLMTS, double Omega, IncField *IF,
               double X[3], cdouble EHTot[6]);

/***************************************************************/
/* displaced surface-integral power, force, torque  ************/
/***************************************************************/
void GetDSIPFT(GLMTSolver *GLMTS, double Omega, IncField *IF,
               double PFT[NUMPFT], double DSIRadius, int DSIPoints);

/***************************************************************/
/* utility routines ********************************************/
/***************************************************************/
cdouble GetEps(GLMTSolver *GLMTS, double Omega, double X[3]);

#endif
