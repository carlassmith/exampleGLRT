/*! \file MakeSubregions.cpp
 *  \author Fang Huang
 *  \date October 10, 2010
 *  \brief This file contains the functions to crop subregions out of the original dataset.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <memory.h>
#include "image_operation.h"
#include "definitions.h"

//#define max(a,b) ( (a) >= (b) ? (a) : (b) )  

void getbox(const float *data,const int ii,const int sz,const int szinX,const int szinY,const int szinZ,const float *X,const float *Y,const float *Z,
			float* dataout,float* leftcoord,float* topcoord);

void MakeSubregions(const int candicc, const float *data, const int szinX, const int szinY, const int szinZ, float *candx,
					float *candy, float *candz, float boxsz, float *dataout, float *left, float *top)
{    
/*! \brief Iterate of the candidates and copy out the subregions.
 *  \param candicc number of candidates
 *  \param data the original dataset
 *  \param szinX data's x dimension
 *  \param szinY data's y dimension
 *  \param szinZ data's z dimension
 *  \param candx array of candidate's x coordinates 
 *  \param candy array of candidate's y coordinates
 *  \param candz array of candidates z corrdinates
 *  \param boxsz size of box to crop out
 *  \param dataout array of copied subregions
 *  \param left location of left edge of the subregion in the original data
 *  \param top  location of the top edge of the subregion in the original data
 */
	int sz,ii;
	
	sz=(int) floor(boxsz); 

	for (ii=0;ii<candicc;ii++)
		getbox(data,ii,sz,(int)szinX,(int)szinY,(int)szinZ,candx,candy,candz,dataout,left,top);

	return;
};

void getbox(const float *data,const int ii,const int sz,const int szinX,const int szinY,const int szinZ,
			const float *X,const float *Y,const float *Z, float* dataout,float* leftcoord,float* topcoord)
{
/*! \brief Copy a subregion of of the original data and return it's left and top location.
 *  \param data the original dataset
 *  \param ii current candidate to copy
 *  \param sz size of the subregion to copy
 *  \param szinX data's x dimension
 *  \param szinY data's y dimension
 *  \param szinZ data's z dimension
 *  \param X array of candidate's x coordinates 
 *  \param Y array of candidate's y coordinates
 *  \param Z array of candidates z corrdinates
 *  \param dataout array of copied subregions
 *  \param leftcoord location of left edge of the subregion in the original data
 *  \param topcoord  location of the top edge of the subregion in the original data
 */
	int yy;
	int l,r,t,b;
	
	const int szl=(int) floor(sz/2.0+.5);
	//get coordinates
	const int x=(int) floor(X[ii]+.5)+1; 
	const int y=(int) floor(Y[ii]+.5)+1;
	const int z=(int) Z[ii];

#ifdef true
	if (x<1 || y<1 || z<0 || x>szinX || y>szinY || z>=szinZ) {
		static char msgbuf[256];
		sprintf(msgbuf,"Point %d out of bounds position %d,%d,%d dataset size %d,%d,%d",ii,x,y,z,szinX, szinY, szinZ);
		mexErrMsgTxt(msgbuf);
	}
#endif

	//mexPrintf("ii: %d %d %d %d %d\n",ii,x,y,z,szl);

	//if (z>szinZ-1)
	//{
	//	mexPrintf("%d %d %d %d\n",x,y,z,szl);
	//    mexErrMsgTxt(" Z index larger than datasize\n");
	//}

	l=max(x-szl,0);
	r=l+sz-1;
	if (r>(szinX-1)) {r=szinX-1;l=r-sz+1;}
	t=max(y-szl,0);
	b=t+sz-1;
	if (b>(szinY-1)) {b=szinY-1;t=b-sz+1;}

	for (yy=0;yy<sz;yy++) {
		//for (xx=0;xx<sz;xx++) dataout[sz*sz*ii+sz*yy+xx]=data[szinX*szinY*z+szinX*(t+yy)+(l+xx)];
		memcpy(dataout+(sz*sz*ii+sz*yy+0),data+(szinX*szinY*z+szinX*(t+yy)+(l+0)),sz*sizeof(float));
	}
	leftcoord[ii]=(float) l;
	topcoord[ii]=(float) t;
	return;
}
