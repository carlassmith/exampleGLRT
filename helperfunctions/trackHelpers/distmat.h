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

	/* Include math.h */
	#include <math.h>

	/* Function prototypes without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. */

	/* calcDistMatrix1D */
	void calcDistMatrix1D(double *D, double *M, double *N, int Mrows, int Nrows);
	
	/* calcDistMatrix2D */
	void calcDistMatrix2D(double *D, double *M, double *N, int Mrows, int Nrows);
	
	/* calcDistMatrix3D */
	void calcDistMatrix3D(double *D, double *M, double *N, int Mrows, int Nrows);

#endif

