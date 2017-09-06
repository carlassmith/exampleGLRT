/*! \file GPUgenerateblobs.cu
 *  \author Fang Huang
 *  \date October 10, 2010
 *  \brief This file contains all relevant code to reassemble the final image
 *  from the fitted subregions.
 */

// includes, system
#include "image_operation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "matrix.h"

// Thread block size
#define BSZ 64
#define MEM 70
#define IMSZ 11
#define IMSZBIG 21
#define imMEM 4000
#define NK 128 //number of blocks to run in each kernel
#define pi 3.141592f
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#define max(a,b)            (((a) > (b)) ? (a) : (b))

void cudasafe( cudaError_t err, char* str, int lineNumber);
void CUDAERROR(const char *instr,int lineNumber);

__global__ void kernel_guassiansampleblobs(int,int,int, float*,float*,float*, float*,float*,float*,float*,float*,float*);
__global__ void kernel_guassianintegrateblobs(int,int,int, float*,float*,float*, float*,float*,float*,float*,float*,float*);

void GPUgenerateblobs(int blobn,int xsz,int ysz,float *xx,float *yy, float *nn, float *LAx, float *LAy, float *cov, float *im)
{
/*!
 *  \brief Setup memory and run the CUDA kernel_guassian methods.
 *  \param blobn 
 *  \param xsz
 *  \param ysz
 *  \param xx 
 *  \param yy
 *  \param nn
 *  \param LAx
 *  \param LAy
 *  \param cov
 *  \param im foo?
 */
	float *d_xarray,*d_yarray,*d_Narray,*d_xsigma,*d_ysigma,*d_covariance,*d_im,*d_xl,*d_yt;
	float *subim,*xl,*yt;
	int ii,iii,jj,kk,loc,boxsz;
	float sigma;
	int BlockSize;
	int flag=0;

	// define text variable
	char text[1024];


	//find max sigma
	float maxsigma=-1;
	for(ii=0;ii<blobn;ii++){
		sigma=sqrt(pow(LAx[ii],2)+pow(LAy[ii],2));
		maxsigma=max(maxsigma,sigma);
	}
	boxsz=(int) round(float (4*maxsigma+1));
	boxsz=min(boxsz,20);// max box size is 20, won't work for psfsigma>6
	//boxsz=20;


	BlockSize=(int) min(ceil((float) 15000/4/boxsz/boxsz),64);
	int memblobsnum=(int)ceil((float)blobn/BlockSize)+128;


	cudaMalloc((void**)&d_xarray, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_xarray, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_xarray, xx, blobn*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_yarray, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_yarray, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_yarray, yy,blobn*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_Narray, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_Narray, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_Narray, nn,blobn*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_xsigma, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_xsigma, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_xsigma, LAx,blobn*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_ysigma, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_ysigma, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_ysigma, LAy,blobn*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_covariance, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_covariance, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_covariance, cov,blobn*sizeof(float), cudaMemcpyHostToDevice);


	cudaMalloc((void**)&d_im, boxsz*boxsz*memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_im, 0, boxsz*boxsz*memblobsnum*BlockSize*sizeof(float));

	cudaMalloc((void**)&d_xl, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_xl, 0, memblobsnum*BlockSize*sizeof(float));

	cudaMalloc((void**)&d_yt, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_yt, 0, memblobsnum*BlockSize*sizeof(float));





	//only run NK blocks in each kernel
	int numK=(int)ceil((float)blobn/BlockSize/NK);

	for (int ii=0;ii<numK;ii++) {

		int blockx = (int) min(ceil(((float)(((float)blobn)/BlockSize)-ii*NK)), NK);
		blockx = max(blockx,1);
		int threadx= BlockSize;


		dim3 dimBlock(threadx);
		dim3 dimGrid(blockx);

		//printf("threadx: %d,blockx: %d\n", threadx, blockx);

		switch (flag)
		{
		case 0:
			kernel_guassiansampleblobs<<<dimGrid, dimBlock>>>(ii,BlockSize,boxsz, d_xarray,d_yarray,d_Narray, d_xsigma,d_ysigma,d_covariance,d_im,d_xl,d_yt);
			CUDAERROR(text,__LINE__);
			break;//15x15 images, 64 per block
		case 1:
			kernel_guassianintegrateblobs<<<dimGrid, dimBlock>>>(ii,BlockSize,boxsz, d_xarray,d_yarray,d_Narray, d_xsigma,d_ysigma,d_covariance,d_im,d_xl,d_yt);
			CUDAERROR(text,__LINE__);
			break;//15x15 images, 64 per block
		}


		// too make the loop work, we have to operate on d_data
		// this trick works to over smart compiler!
		cudaMemcpy(d_xarray, xx, 1*sizeof(float), cudaMemcpyHostToDevice);

		//printf("Fitting subimages: %d/%d is DONE...\n", (ii+1), numK);
		//CUDAERRROR("kernel");
		mexEvalString("pause(0.001)");

	}

	subim= (float * )malloc(blobn*boxsz*boxsz*sizeof(float));
	memset (subim,0,blobn*boxsz*boxsz*sizeof(float));
	xl=(float * )malloc(blobn*sizeof(float));
	memset(xl,0,blobn*sizeof(float));
	yt=(float * )malloc(blobn*sizeof(float));
	memset(yt,0,blobn*sizeof(float));

	//reconstruct images

	cudaMemcpy(subim, d_im, blobn*boxsz*boxsz*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(xl, d_xl, blobn*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(yt, d_yt, blobn*sizeof(float), cudaMemcpyDeviceToHost);

	for(kk=0;kk<blobn;kk++){
		for(jj=0;jj<boxsz;jj++){
			for(iii=0;iii<boxsz;iii++){
				if ((((int)xl[kk]+iii)<xsz)&&(((int)yt[kk]+jj)<ysz)){
					loc=((int)yt[kk]+jj)*xsz+(int)xl[kk]+iii;
					if((subim[kk*boxsz*boxsz+jj*boxsz+iii]>0)&&(subim[kk*boxsz*boxsz+jj*boxsz+iii]<100000))
					im[loc]+=subim[kk*boxsz*boxsz+jj*boxsz+iii];	
				}
			}
		}
	}




	free(subim);
	free(xl);
	free(yt);
	cudaFree(d_xarray);
	cudaFree(d_yarray);
	cudaFree(d_Narray);
	cudaFree(d_xsigma);
	cudaFree(d_ysigma);
	cudaFree(d_covariance);
	cudaFree(d_im);
	cudaFree(d_xl);
	cudaFree(d_yt);


}


//kernel_guassiansampleblobs<<<dimGrid, dimBlock>>>(ii,blockx,BlockSize,sz, d_xarray,d_yarray,d_Narray, d_xsigma,d_ysigma,d_covariance,d_im,d_xl,d_yt);   //15x15 images, 64 per block

__global__ void kernel_guassiansampleblobs(int iiK,int BlockSize, int sz, float *d_xarray,float *d_yarray,
										   float *d_Narray, float *d_xsigma,float *d_ysigma,float *d_covariance,
										   float *d_im,float *d_xl,float *d_yt  ) 
{
/*!
 *  \brief dunno.
 *  \param iiK
 *  \param BlockSize
 *  \param sz
 *  \param d_xarray 
 *  \param d_yarray
 *  \param d_Narray
 *  \param d_xsigma
 *  \param d_ysigma
 *  \param d_covariance
 *  \param d_im
 *  \param d_xl
 *  \param d_yt
 */
	int tx = threadIdx.x; //matrix number index
	int bx = blockIdx.x;
	float x,y,xsigma,ysigma,covariance,N;
	int xl;
	int yt;
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
	covariance=d_covariance[bx*BlockSize+tx]/xsigma/ysigma;

	xl=round(x)-round(float(sz/2-1) );
	xl=max(xl,0);

	yt=round(y)-round(float (sz/2-1));
	yt=max(yt,0);


	for (ii=0;ii<sz;ii++) 
	{
		for(jj=0;jj<sz;jj++) {

			// generate model for pixel ii jj
			pixelx=ii;
			pixely=jj;
			s_im[tx*sz*sz+jj*sz+ii]=N/(2*pi*xsigma*ysigma*sqrt(1-pow(covariance,2)))*exp(-1/(2*(1-pow(covariance,2)))*(pow(x-xl-pixelx,2)/pow(xsigma,2)+pow(y-yt-pixely,2)/pow(ysigma,2)-2*covariance*(x-xl-pixelx)*(y-yt-pixely)/(xsigma*ysigma)));
		}
	}



	for (ii=0;ii<sz;ii++){
		for(jj=0;jj<sz;jj++)
		{
			d_im[bx*BlockSize*sz*sz+tx*sz*sz+jj*sz+ii]=s_im[tx*sz*sz+jj*sz+ii];
			d_xl[bx*BlockSize+tx]=xl;
			d_yt[bx*BlockSize+tx]=yt;
		}
	}
	return;



}



__global__ void kernel_guassianintegrateblobs(int iiK,int BlockSize, int sz, float *d_xarray,float *d_yarray,
											  float *d_Narray, float *d_xsigma,float *d_ysigma,float *d_covariance,
											  float *d_im,float *d_xl,float *d_yt  ) 
{
/*!
 *  \brief dunno.
 *  \param iiK
 *  \param BlockSize
 *  \param sz
 *  \param d_xarray 
 *  \param d_yarray
 *  \param d_Narray
 *  \param d_xsigma
 *  \param d_ysigma
 *  \param d_covariance
 *  \param d_im
 *  \param d_xl
 *  \param d_yt
 */
	int tx = threadIdx.x; //matrix number index
	int bx = blockIdx.x;
	float x,y,xsigma,ysigma,covariance;
	//float covariance
	float N;
	int xl;
	int yt;
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
	covariance=d_covariance[bx*BlockSize+tx];

	xl=round(x)-round(float (sz/2-1));
	xl=max(xl,0);

	yt=round(y)-round(float (sz/2-1));
	yt=max(yt,0);


	for (ii=0;ii<sz;ii++) {
		for(jj=0;jj<sz;jj++) {

			// generate model for pixel ii jj
			pixelx=ii;
			pixely=jj;
			s_im[tx*sz*sz+jj*sz+ii]=N/4*(erf((x-xl-pixelx-0.5)/sqrt(2*pow(xsigma,2)))-erf((x-xl-pixelx+0.5)/sqrt(2*pow(xsigma,2))))*(erf((y-yt-pixely-0.5)/sqrt(2*pow(ysigma,2)))-erf((y-yt-pixely+0.5)/sqrt(2*pow(ysigma,2))));  //exp(-1/(2*(1-pow(covariance,2)))*(pow(x-xl-pixelx,2)/pow(xsigma,2)+pow(y-yt-pixely,2)/pow(ysigma,2)-2*covariance*(x-xl-pixelx)*(y-yt-pixely)/(xsigma*ysigma)));
		}
	}


	for (ii=0;ii<sz;ii++) {
		for(jj=0;jj<sz;jj++)
		{
			d_im[bx*BlockSize*sz*sz+tx*sz*sz+jj*sz+ii]=s_im[tx*sz*sz+jj*sz+ii];
			d_xl[bx*BlockSize+tx]=xl;
			d_yt[bx*BlockSize+tx]=yt;
		}
	}
	return;



}








