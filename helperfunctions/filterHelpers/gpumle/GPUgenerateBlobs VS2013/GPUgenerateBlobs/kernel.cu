
#include "cuda_runtime.h"
#include "definitions.h"
#include "kernel.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

__global__ void kernel_guassiansampleblobs(int iiK,int BlockSize, int sz, float *d_xarray,float *d_yarray,float *d_Narray, float *d_xsigma,float *d_ysigma,float *d_covariance,float *d_im,float *d_xl,float *d_yt  ) {
	int tx = threadIdx.x; //matrix number index
	int bx = blockIdx.x;
	float x,y,xsigma,ysigma,covariance,N;
	float xl;
	float yt;
	int ii,jj,pixelx,pixely;


//	float model;//

	__shared__ float s_im[imMEM];


	bx=bx+iiK*NK;
	//import datas from device to shared memory

	x=d_xarray[bx*BlockSize+tx];
	y=d_yarray[bx*BlockSize+tx];
	N=d_Narray[bx*BlockSize+tx];
	xsigma=d_xsigma[bx*BlockSize+tx];
	ysigma=d_ysigma[bx*BlockSize+tx];
	covariance=d_covariance[bx*BlockSize+tx];
	xl=round(x)-round(float (sz/2-1));
	xl=max(xl,0.0);

	yt=round(y)-round(float (sz/2-1));
	yt=max(yt,0.0);


	for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {

		// generate model for pixel ii jj
		pixelx=ii;
		pixely=jj;
		s_im[tx*sz*sz+jj*sz+ii]=N/(2*M_PI*xsigma*ysigma*sqrt(1-pow(covariance,2)))*exp(-1/(2*(1-pow(covariance,2)))*(pow(x-xl-pixelx,2)/pow(xsigma,2)+pow(y-yt-pixely,2)/pow(ysigma,2)-2*covariance*(x-xl-pixelx)*(y-yt-pixely)/(xsigma*ysigma)));
	}



	for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++)
	{
		d_im[bx*BlockSize*sz*sz+tx*sz*sz+jj*sz+ii]=s_im[tx*sz*sz+jj*sz+ii];
		d_xl[bx*BlockSize+tx]=xl;
		d_yt[bx*BlockSize+tx]=yt;
	}

// 	return true;



}



__global__ void kernel_guassianintegrateblobs(int iiK,int BlockSize, int sz, float *d_xarray,float *d_yarray,float *d_Narray, float *d_xsigma,float *d_ysigma,float *d_covariance,float *d_im,float *d_xl,float *d_yt  ) {
	int tx = threadIdx.x; //matrix number index
	int bx = blockIdx.x;
	float x,y,xsigma,ysigma,N;
	float xl;
	float yt;
	int ii,jj,pixelx,pixely;


	//float model;//

	__shared__ float s_im[imMEM];


	bx=bx+iiK*NK;
	//import datas from device to shared memory

	x=d_xarray[bx*BlockSize+tx];
	y=d_yarray[bx*BlockSize+tx];
	N=d_Narray[bx*BlockSize+tx];
	xsigma=d_xsigma[bx*BlockSize+tx];
	ysigma=d_ysigma[bx*BlockSize+tx];
	//covariance=d_covariance[bx*BlockSize+tx];

	xl=round(x)-round(float (sz/2-1));
	xl=max(xl,0.0);

	yt=round(y)-round(float (sz/2-1));
	yt=max(yt,0.0);

	for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {

		// generate model for pixel ii jj
		pixelx=ii;
		pixely=jj;
		s_im[tx*sz*sz+jj*sz+ii]=N/4*(erf((x-xl-pixelx-0.5)/sqrt(2*pow(xsigma,2)))-erf((x-xl-pixelx+0.5)/sqrt(2*pow(xsigma,2))))*(erf((y-yt-pixely-0.5)/sqrt(2*pow(ysigma,2)))-erf((y-yt-pixely+0.5)/sqrt(2*pow(ysigma,2)))); 
	}



	for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++)
	{
		d_im[bx*BlockSize*sz*sz+tx*sz*sz+jj*sz+ii]=s_im[tx*sz*sz+jj*sz+ii];
		d_xl[bx*BlockSize+tx]=xl;
		d_yt[bx*BlockSize+tx]=yt;
	}

// 	return;



}