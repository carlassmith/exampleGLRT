/* -------------------------------------------------------    */
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

#ifndef DISTMAT_H_
    
    #define DISTMAT_H_

    #include "distmat.h"

#endif
        
#include <math.h>        

/* calcDistMatrix1D */

void calcDistMatrix1D(double *D, double *M, double *N, int Mrows, int Nrows)
{
    int i,j;
    
    for (i=0;i<Nrows;i++) {
        
        for (j=0;j<Mrows;j++) {
        
            // Store calculated distance in D
            *D++=*(N+i)-*(M+j);            
        }
    }
}

/* calcDistMatrix2D */

void calcDistMatrix2D(double *D, double *M, double *N, int Mrows, int Nrows)
{
    double	mX=0, mY=0;
    double	nX=0, nY=0;
    double  *posM, *posN;
    int  i,j;
    
    for (i=0;i<Nrows;i++) {
        
        // Get source position
        posN=N+i;
        nX=*posN;
        nY=*(posN+Nrows);
        
        for (j=0;j<Mrows;j++) {
        
            // Get target position
            posM=M+j;
            mX=*posM;
			mY=*(posM+Mrows);
            
             // Store calculated distance in D
            *D++=sqrt((nY-mY)*(nY-mY)+(nX-mX)*(nX-mX));           
        }
    }
}

/* calcDistMatrix3D */

void calcDistMatrix3D(double *D, double *M, double *N, int Mrows, int Nrows)
{
    double	mX=0, mY=0, mZ=0;
    double	nX=0, nY=0, nZ=0;
    double  *posM, *posN;
    int  i,j;
    
    for (i=0;i<Nrows;i++) {
        
        // Get source position
        posN=N+i;
        nX=*posN;
        nY=*(posN+Nrows);
		nZ=*(posN+2*Nrows);
        
        for (j=0;j<Mrows;j++) {
        
            // Get target position
            posM=M+j;
            mX=*posM;
			mY=*(posM+Mrows);
			mZ=*(posM+2*Mrows);
            
            // Store calculated distance in D
            *D++=sqrt((nY-mY)*(nY-mY)+(nX-mX)*(nX-mX)+(nZ-mZ)*(nZ-mZ));
            
        }
    }
}
