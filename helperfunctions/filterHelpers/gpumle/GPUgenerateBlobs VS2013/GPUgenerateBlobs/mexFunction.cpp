
#include <stdio.h> /* for stderr */
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include "mex.h"
#include "math.h"

#include "cuda_runtime.h"

#include "matrix.h"
#include "definitions.h"
#include "kernel.h"

extern void kernel_guassiansampleblobs_wrapper(dim3 dimGrid, dim3 dimBlock, int,int,int, float*,float*,float*, float*,float*,float*,float*,float*,float*);
extern void kernel_guassianintegrateblobs_wrapper(dim3 dimGrid, dim3 dimBlock,int,int,int, float*,float*,float*, float*,float*,float*,float*,float*,float*);

void CUDAERROR(const char *instr);
void cudasafe( cudaError_t err, char* str, int lineNumber);

#define min(a,b)            (((a) < (b)) ? (a) : (b))
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#define roundf(x)           x >= 0.0f ? floorf(x + 0.5f) : ceilf(x - 0.5f);


//*******************************************************************************************
void CUDAERROR(const char *instr) {
	cudaError_t errornum;
	const char *str;
	if (errornum=cudaGetLastError()) {
		str=cudaGetErrorString(errornum);
		mexPrintf("%s\n", str);
		mexPrintf("You should clear this function in MATLAB for proper operation.\n", str);
	}
}


//*******************************************************************************************
void mexFunction(int nlhs, mxArray *plhs[],	int	nrhs, const	mxArray	*prhs[]) {
/*!
 *  \brief Entry point in the code for Matlab.  Equivalent to main().
 *  \param nlhs number of left hand mxArrays to return
 *  \param plhs array of pointers to the output mxArrays
 *  \param nrhs number of input mxArrays
 *  \param prhs array of pointers to the input mxArrays.
 */
	int blockx;
	int threadx;
	int ii,iii,jj,kk,flag;
	int memblobsnum,ysz,xsz;
	float * xarray, * yarray, * Narray,*yt,*xl,*xsigma,*ysigma,*covariance,*im;
	float *d_xarray, *d_yarray, *d_Narray, *d_xsigma, *d_ysigma,*d_covariance,*d_im,*d_xl,*d_yt,*subim;
	const mwSize *datasize;
	int locr;
	mwSize imdim[2];



	if (nrhs<9)
		mexErrMsgTxt("xsize,ysize,x_array, y_array, N_array, sigmaX, sigmaY, covariance, UseIntegrated_FLAG\n");

	if (mxGetClassID(prhs[0])!=mxSINGLE_CLASS)
		mexErrMsgTxt("Data must be comprised of single floats!\n");
	if (mxGetClassID(prhs[1])!=mxSINGLE_CLASS)
		mexErrMsgTxt("Data must be comprised of single floats!\n");
	if (mxGetClassID(prhs[2])!=mxSINGLE_CLASS)
		mexErrMsgTxt("Data must be comprised of single floats!\n");
	if (mxGetClassID(prhs[3])!=mxSINGLE_CLASS)
		mexErrMsgTxt("Data must be comprised of single floats!\n");
	if (mxGetClassID(prhs[4])!=mxSINGLE_CLASS)
		mexErrMsgTxt("Data must be comprised of single floats!\n");
	if (mxGetClassID(prhs[5])!=mxSINGLE_CLASS)
		mexErrMsgTxt("Data must be comprised of single floats!\n");


	datasize=mxGetDimensions(prhs[2]);
	if (datasize[1]!=1)
		mexErrMsgTxt("xarray should be n X 1 array\n");

	datasize=mxGetDimensions(prhs[3]);

	if (datasize[1]!=1)
		mexErrMsgTxt("xarray should be n X 1 array\n");

	datasize=mxGetDimensions(prhs[4]);

	if (datasize[1]!=1)
		mexErrMsgTxt("xarray should be n X 1 array\n");

	datasize=mxGetDimensions(prhs[5]);
	if (datasize[1]!=1)
		mexErrMsgTxt("xarray should be n X 1 array\n");


	xsz =(int) mxGetScalar(prhs[0]);
	ysz =(int) mxGetScalar(prhs[1]);
	imdim[0]=xsz;
	imdim[1]=ysz;
	//PSFSigma=(float)mxGetScalar(prhs[1]); //matlab-dip_image convention
	xarray =(float *) mxGetData(prhs[2]);
	yarray =(float *) mxGetData(prhs[3]);
	Narray =(float *) mxGetData(prhs[4]);
	xsigma =(float *)mxGetData(prhs[5]);
	ysigma =(float *)mxGetData(prhs[6]);
	covariance =(float *)mxGetData(prhs[7]);
	flag =(int) mxGetScalar(prhs[8]);


	int blobn = (int)datasize[0];
	//check to see whether it exceeds 400MB=20*20*250000*4, if it reaches, it probably could crush. The image reconstruction should be implemented on reconstruct on GPU device memory instead of host memory on PC.
	if (blobn>250000)
		mexErrMsgTxt("Size of total blobs to be generated exceeded capacity of GPU global memory (blobn>250,000). Consider partitioning them using 'for' loops \n");

	float maxsigma=-1;
	float sigma;
	for(ii=0;ii<blobn;ii++){
		sigma=sqrt(pow(xsigma[ii],2)+pow(ysigma[ii],2));
		maxsigma=max(maxsigma,sigma);
	}

	int sz= (int) roundf(8.0*maxsigma);
	sz=min(sz,20);


	if ((flag!=1)&&(flag!=0))
		mexErrMsgTxt("flag can only be 0 or 1\n");

	// over allocate for additional thread reading error
	int BlockSize=(int) min(ceil((float) 15000/4/sz/sz),64);
	memblobsnum=(int)ceil((float)datasize[0]/BlockSize)+128;


	cudaMalloc((void**)&d_xarray, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_xarray, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_xarray, xarray, datasize[0]*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_yarray, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_yarray, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_yarray, yarray,datasize[0]*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_Narray, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_Narray, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_Narray, Narray,datasize[0]*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_xsigma, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_xsigma, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_xsigma, xsigma,datasize[0]*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_ysigma, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_ysigma, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_ysigma, ysigma,datasize[0]*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_covariance, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_covariance, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_covariance, covariance,datasize[0]*sizeof(float), cudaMemcpyHostToDevice);


	cudaMalloc((void**)&d_im, sz*sz*memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_im, 0, sz*sz*memblobsnum*BlockSize*sizeof(float));

	cudaMalloc((void**)&d_xl, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_xl, 0, memblobsnum*BlockSize*sizeof(float));

	cudaMalloc((void**)&d_yt, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_yt, 0, memblobsnum*BlockSize*sizeof(float));


	//only run NK blocks in each kernel
	int numK=(int)ceil((float)datasize[0]/BlockSize/NK);

	for (int ii=0;ii<numK;ii++) {

		blockx = (int) min(ceil(((float)(((float)datasize[0])/BlockSize)-ii*NK)), NK);
		blockx = (int) max(blockx,1);
		threadx= BlockSize;


		dim3 dimBlock(threadx);
		dim3 dimGrid(blockx);

		//printf("threadx: %d,blockx: %d\n", threadx, blockx);

		switch (flag)
		{
		case 0:
			kernel_guassiansampleblobs_wrapper(dimGrid, dimBlock,ii,BlockSize,sz, d_xarray,d_yarray,d_Narray, d_xsigma,d_ysigma,d_covariance,d_im,d_xl,d_yt);
			break;//15x15 images, 64 per block
		case 1:
			kernel_guassianintegrateblobs_wrapper(dimGrid, dimBlock,ii,BlockSize,sz, d_xarray,d_yarray,d_Narray, d_xsigma,d_ysigma,d_covariance,d_im,d_xl,d_yt);
			break;//15x15 images, 64 per block
		}


		// too make the loop work, we have to operate on d_data
		// this trick works to over smart compiler!
		//cudaMemcpy(d_xarray, xarray, 1*sizeof(float), cudaMemcpyHostToDevice);


		CUDAERROR("kernel");
		//mexEvalString("pause(0.001)");

	}

	subim= (float * )malloc(datasize[0]*sz*sz*sizeof(float));
	xl=(float * )malloc(datasize[0]*sizeof(float));
	yt=(float * )malloc(datasize[0]*sizeof(float));


	//reconstruct images
	plhs[0]=mxCreateNumericArray(2, imdim, mxSINGLE_CLASS, mxREAL);
	im=(float *)mxGetData(plhs[0]);

	cudaMemcpy(subim, d_im, datasize[0]*sz*sz*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(xl, d_xl, datasize[0]*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(yt, d_yt, datasize[0]*sizeof(float), cudaMemcpyDeviceToHost);


	for(kk=0;kk<blobn;kk++){
		for(jj=0;jj<sz;jj++){
			for(iii=0;iii<sz;iii++){
				if ((((int)xl[kk]+iii)<(xsz))&&(((int)yt[kk]+jj)<(ysz))){
					locr=((int)yt[kk]+jj)*xsz+(int)xl[kk]+iii;
					if((subim[kk*sz*sz+jj*sz+iii]>0)&&(subim[kk*sz*sz+jj*sz+iii]<100000)&&(locr>=0)&&(locr<((xsz)*(ysz))))
					im[locr]+=subim[kk*sz*sz+jj*sz+iii];	
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