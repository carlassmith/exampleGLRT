/* -------------------------------------------------------    */
/*                                                            */
/* createDistanceMatrix [MATLAB C-MEX]                        */
/*                                                            */
/* -------------------------------------------------------    */
/*                                                            */
/* Files:                                                     */
/*                                                            */
/*     createDistanceMatrix.c - MEX interface                 */
/*     distmat.h              - function prototypes           */
/*     distmat.c              - function definitions          */
/*                                                            */
//* -------------------------------------------------------   */
/*                                                            */
/* distmat.h [distance matrix - function prototypes]          */
/*                                                            */
/*                                                            */ 
/* This file is part of u-track.                              */
/*                                                            */ 
/*    u-track is free software: you can redistribute it       */ 
/*    and/or modify it under the terms of the GNU General     */
/*    Public License as published by the Free Software        */ 
/*    Foundation, either version 3 of the License, or         */
/*    (at your option) any later version.                     */  
/*                                                            */
/*   u-track is distributed in the hope that it will be       */ 
/*   useful, but WITHOUT ANY WARRANTY; without even           */ 
/*   the implied warranty of without even the implied         */ 
/*   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR  */ 
/*   PURPOSE.GNU General Public License for more details.     */
/*                                                            */
/*    You should have received a copy of the GNU General      */
/*    Public License along with u-track.                      */
/*    If not, see <http://www.gnu.org/licenses/>.             */
/*                                                            */
/*                                                            */ 
/*                                                            */
/* Copyright: Aaron Ponti - 02/08/28                          */
/*                                                            */
/* -------------------------------------------------------    */


#include "mex.h"
#include "distmat.h"

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    /* Initialize pointers for 2 inputs and 1 output */
    double *M, *N;
	double *D;

	/* Initialize int variables to store matrix dimensions */
	int Mrows,Mcols;
	int Nrows,Ncols;

    /* Check that the number of input and output parameters is valid */
	if(nrhs != 2)
		mexErrMsgTxt("Two input parameters required.");
	if(nlhs > 1)
		mexErrMsgTxt("One output parameter required.");	
	
	/* Read input parameter dimensions */
	Mrows=mxGetM(prhs[0]);
	Mcols=mxGetN(prhs[0]);
	Nrows=mxGetM(prhs[1]);
	Ncols=mxGetN(prhs[1]);
	
	/* Check input parameter dimension */
	if ((Mcols>3) || (Ncols>3))
		mexErrMsgTxt("Point coordinates in more than 3 dimensions are not supported.");
	
	if (Mcols!=Ncols)
		mexErrMsgTxt("The points in the coordinate matrices have different number of dimensions.");

	/* Create matrix for the return argument D */
	plhs[0]=mxCreateDoubleMatrix(Mrows,Nrows, mxREAL);
    
    /* Assign pointers to each input and output */
	M=mxGetPr(prhs[0]);
	N=mxGetPr(prhs[1]);
	D=mxGetPr(plhs[0]);

   	/* Call the correcponding C function */
	if (Mcols==1) { calcDistMatrix1D(D,M,N,Mrows,Nrows); }
	if (Mcols==2) { calcDistMatrix2D(D,M,N,Mrows,Nrows); }
	if (Mcols==3) { calcDistMatrix3D(D,M,N,Mrows,Nrows); }
	
}
