#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "definitions.h"
#include "image_operation.h"
#include "GPUgaussLib.cuh"
#define BSZ 64
#define MEM 3872
#define IMSZ 11
#define IMSZBIG 21
#define NK 128 //number of blocks to run in each kernel
#define pi 3.141592f
#define NUM 8
#define a1 -0.9911f
#define a2 -0.6763f
#define b1 0.8055f
#define b2 -1.2451f

__global__ void kernel_MLEFit(int iiK, float *d_data, float PSFSigma, int sz, int iterations, int BlockSize,
							  int num, float Nave, float Nsigma, float llThreshold, float *d_X, float *d_Y, 
							  float *d_N, float *d_b, float *d_Div, float * d_Dett);

void cudasafe( cudaError_t err, char* str, int lineNumber);
void CUDAERROR(const char *instr,int lineNumber);

void GPUmultifit(int Nfitraw, float *in_sub,const int sz, float psfsigma,const int iterations,int fitnum, float Nave, 
				 const float Nsigmafit,float threshold,float *xx,float *yy,float *nn,float *bb,float *div)
{
	/*!
	*  \brief Setup and run CUDA kernel_MLEFit.
	*  \param NFitraw
	*  \param in_sub
	*  \param sz
	*  \param psfsigma 
	*  \param iterations
	*  \param fitnum
	*  \param Nave
	*  \param Nsigmafit
	*  \param threshold
	*  \param xx
	*  \param yy
	*  \param nn
	*  \param bb
	*  \param div
	*/
	int blockx;
	int threadx;
	//int iiK,
	int Nfits;
	float *d_data,*d_X, *d_Y,* d_N,* d_b,*d_Div,* d_Dett;


	int BlockSize = (int) floor((float)15000/4/sz/sz);
	BlockSize = max(4, BlockSize);
	BlockSize = min(BSZ, BlockSize);
	Nfits=(int) floor(BlockSize*ceil( (float) Nfitraw/BlockSize));

	cudaMalloc((void**)&d_data, sz*sz*Nfits*sizeof(float));
	cudaMemset(d_data, 0, sz*sz*Nfits*sizeof(float));
	cudaMemcpy(d_data, in_sub, sz*sz*Nfitraw*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_X,    fitnum*Nfits*sizeof(float));
	cudaMemset(d_X, 0, fitnum*Nfits*sizeof(float));
	cudaMalloc((void**)&d_Y,    fitnum*Nfits*sizeof(float));
	cudaMemset(d_Y, 0, fitnum*Nfits*sizeof(float));
	cudaMalloc((void**)&d_b,    Nfits*sizeof(float));
	cudaMemset(d_b,0,Nfits*sizeof(float));
	cudaMalloc((void**)&d_N,    fitnum*Nfits*sizeof(float));
	cudaMemset(d_N, 0, fitnum*Nfits*sizeof(float));
	cudaMalloc((void**)&d_Div,  Nfits*sizeof(float));
	cudaMemset(d_Div, 0, Nfits*sizeof(float));
	cudaMalloc((void**)&d_Dett, Nfits*sizeof(float));
	cudaMemset(d_Dett, 0, Nfits*sizeof(float));

	int numK=(int) ceil((float)Nfits/BlockSize/NK);

	// declare text variable
	char text[1024];

	printf("Fitting subimages\n");

	for (int ii=0;ii<numK;ii++) 
	{
		//setup kernel
		blockx = min(Nfits/BlockSize-ii*NK, NK);
		threadx= BlockSize;

		dim3 dimBlock(threadx);
		dim3 dimGrid(blockx);

		printf("threadx: %d,blockx: %d,Nfitraw: %d\n", threadx, blockx, Nfitraw);

		kernel_MLEFit<<<dimGrid, dimBlock>>>(ii, d_data, psfsigma, sz, iterations, BlockSize, fitnum, Nave, Nsigmafit, threshold, d_X, d_Y, d_N, d_b,d_Div, d_Dett);
		CUDAERROR(text,__LINE__);

		// too make the loop work, we have to operate on d_data
		// this trick works to over smart compiler!

		cudaMemcpy(d_data, in_sub, 1*sizeof(float), cudaMemcpyHostToDevice);
		printf("Fitting subimages: %d/%d is DONE...\n", (ii+1), numK);

		mexEvalString("pause(0.001)");



	}

	cudaMemcpy(xx, d_X, Nfitraw*fitnum*sizeof(float), cudaMemcpyDeviceToHost); //reversed for dipimage

	cudaMemcpy(yy, d_Y, Nfitraw*fitnum*sizeof(float), cudaMemcpyDeviceToHost);//reversed for dipimage

	cudaMemcpy(nn, d_N, Nfitraw*fitnum*sizeof(float), cudaMemcpyDeviceToHost);

	cudaMemcpy(bb, d_b, Nfitraw*sizeof(float), cudaMemcpyDeviceToHost);

	cudaMemcpy(div, d_Div, Nfitraw*sizeof(float), cudaMemcpyDeviceToHost);

	//
	//cleanup
	cudaFree(d_X);
	cudaFree(d_Y);
	cudaFree(d_N);
	cudaFree(d_b);
	cudaFree(d_Div);
	cudaFree(d_data);
	cudaFree(d_Dett);


}

__global__ void kernel_MLEFit(int iiK, float *d_data, float PSFSigma, int sz, int iterations, int BlockSize,
							  int num, float Nave, float Nsigma, float llThreshold, float *d_X, float *d_Y, 
							  float *d_N, float *d_b, float *d_Div, float * d_Dett) 
{
	/*!
	*  \brief Gaussian fit the pixels location in the various subregions.
	*  \param iiK
	*  \param d_data
	*  \param PSFSigma
	*  \param sz 
	*  \param iterations
	*  \param BlockSize
	*  \param num
	*  \param Nave
	*  \param Nsigma
	*  \param llThreshold
	*  \param d_X
	*  \param d_Y
	*  \param d_N
	*  \param d_b
	*  \param d_Div
	*  \param d_Dett
	*/
	__shared__ float s_data[MEM];

	float Nmax;

	int tx = threadIdx.x; //fits
	int bx = blockIdx.x;
	int ii, jj, kk, qq, mm, nn, nnM, goodness;
	int nnfit=0;
	//int nummax;
	int breaktmp, nncount, rimsigny, rimsignx;
	float dLLx, dLLy, dLLN, dLLb;
	float deflamodel;
	float tmp, signx, signy, meanx, meany;
	float b, N, x, y, xmin, xmax, ymin, ymax;
	float model, cf, df, data,bini;
	float norm;
	float Div;
	float bfit=0;
	//float Dett;
	//float Ixx, Iyy, III, Ibgbg;
	//float Ixy, IyI, Iybg, IIbg;
	//float Ixbg, IxI;
	float imdy, imdx, imdI, /*imdbg,*/ imddx, imddy, imddI;
	float PSFy, PSFx, pval, zval,minDiv,maxpval;
	float stepxtot;
	float stepytot;
	float stepbgtot;
	float stepItot;
	float tmpx, tmpy, tmpsum;
	float xarray[NUM]={0, 0, 0, 0, 0, 0, 0, 0};
	float yarray[NUM]={0, 0, 0, 0, 0, 0, 0, 0};
	float Narray[NUM]={0, 0, 0, 0, 0, 0, 0, 0};
	float xarrayfit[NUM]={0, 0, 0, 0, 0, 0, 0, 0}, yarrayfit[NUM]={0, 0, 0, 0, 0, 0, 0, 0}, Narrayfit[NUM]={0, 0, 0, 0, 0, 0, 0, 0};
	float x_current[NUM]={0, 0, 0, 0, 0, 0, 0, 0}, y_current[NUM]={0, 0, 0, 0, 0, 0, 0, 0}, N_current[NUM]={0, 0, 0, 0, 0, 0, 0, 0};


	//increase index because of looping of kernels Preventing crash?
	bx=bx+iiK*NK;

	//load data
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++)
        s_data[sz*sz*tx+sz*jj+ii]=d_data[sz*sz*bx*BlockSize+sz*sz*tx+sz*jj+ii];
	
	//initial values
	kernel_CenterofMass2D(sz, &s_data[sz*sz*tx], &x, &y);
	kernel_GaussFMaxMin2D(sz, PSFSigma, &s_data[sz*sz*tx], &Nmax, &b);
	bini=max(b, 0.0001f);

	N = max(0.0f, (Nmax-b)*2*pi*PSFSigma*PSFSigma);

	//N=(float)maxN[tx]*2.0f*pi*PSFSigma*PSFSigma;
	PSFy = 0.0f; PSFx=0.0f;
	//S=(float)PSFSigma;
	norm=1.0f/2.0f/PSFSigma/PSFSigma;
	//start multifit from single fit til num-fit
	tmp=Nave/2/pi/2/PSFSigma/PSFSigma;
	maxpval=-1.0f;
	Nsigma=0.00001;//fixed intensity for first round iteration loop
	for (mm=1;mm<=num;mm++) {
		breaktmp=0;
		goodness=0;
		b=bini;

		//INITIAL GUESSES

		//if fitting for a single fluorophore, use center of mass as initial guess
		if (mm==1){
			xarray[0]=x;
			yarray[0]=y;
			Narray[0]=Nave;

		}
		//if fitting more than one fluorophore, use deflation method and find max
		else {
			tmp=Nave/2/pi/2/PSFSigma/PSFSigma;

			for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
				model=b;
				nncount=0;
				meanx=0;
				meany=0;
				for (nn=0;nn<mm-1;nn++) {
					x=xarray[nn];
					y=yarray[nn];
					N=Narray[nn];
					meanx+=xarray[nn];
					meany+=yarray[nn];
					nncount=nncount+1;
					PSFx=kernel_IntGauss1D(ii, x, PSFSigma);
					PSFy=kernel_IntGauss1D(jj, y, PSFSigma);					
					model+=N*PSFx*PSFy;
				}
				deflamodel=s_data[sz*sz*tx+sz*jj+ii]-model;
				signx=0;
				signy=0;
				rimsignx=1;
				rimsigny=1;
				if (deflamodel>tmp) {
					if (ii==0) rimsignx=-1;
					if (ii==sz-1) rimsignx=-1;
					if (jj==0) rimsigny=-1;
					if (jj==sz-1) rimsigny=-1;

					signx=((ii-meanx/nncount) > 0)? 1 : -1;
					signy=((jj-meany/nncount) > 0)? 1 : -1;
					xarray[mm-1]=ii-0.5f*signx*rimsignx;
					yarray[mm-1]=jj-0.5f*signy*rimsigny;
					Narray[mm-1]=Nave;
					tmp=deflamodel;

				}
			}
			if (tmp==Nave/2/pi/2/PSFSigma/PSFSigma)
				breaktmp=1;

		}

		if (breaktmp==1)
			break;


		//MAIN ITERATIVE LOOP

		for (kk=0;kk<iterations;kk++){

			dLLb=0.0f;
			stepbgtot = 0;

			//generating the model and calc b update
			for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {

				model=b;
				for (nn=0;nn<mm;nn++) {
					x_current[nn]=xarray[nn];
					y_current[nn]=yarray[nn];
					N_current[nn]=Narray[nn];

					x=xarray[nn];
					y=yarray[nn];
					N=Narray[nn];

					PSFx=kernel_IntGauss1D(ii, x, PSFSigma);
					PSFy=kernel_IntGauss1D(jj, y, PSFSigma);
					model+=N*PSFx*PSFy;
				}
				data=s_data[sz*sz*tx+sz*jj+ii];
				cf=0.0f;
				df=0.0f;
				if (model>10e-3f) cf=data/model-1;
				if (model>10e-3f) df=data/pow(model, 2);
				cf=min(cf, 10e4f);
				df=min(df, 10e4f);
				stepbgtot += -df;
				dLLb+=cf;
			}
			b-=min(max(dLLb/stepbgtot, -1e0f), 1e0f)/mm/2;
			b=max(b, 0.001f);
			//b=0.01f;



			//This starts iterative routine for theta_i other than theta_bg which is calculated above.
			for (nn=0;nn<mm;nn++){


				dLLx=0.0f;
				dLLy=0.0f;
				dLLN=0.0f;
				stepxtot = 0.0f;
				stepytot = 0.0f;
				stepItot = 0.0f;

				for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
					//generate model using initial value or iteration value.
					data=s_data[sz*sz*tx+sz*jj+ii];
					model=b;
					for (nnM=0;nnM<mm;nnM++) {
						x=x_current[nnM];
						y=y_current[nnM];
						N=N_current[nnM];
						PSFx=kernel_IntGauss1D(ii, x, PSFSigma);
						PSFy=kernel_IntGauss1D(jj, y, PSFSigma);
						model+=N*PSFx*PSFy;
					}

					cf=0.0f;
					df=0.0f;
					if (model>10e-3f) cf=data/model-1;
					if (model>10e-3f) df=data/pow(model, 2);
					cf=min(cf, 10e4f);
					df=min(df, 10e4f);

					x=x_current[nn];
					y=y_current[nn];
					N=N_current[nn];

					PSFx=kernel_IntGauss1D(ii, x, PSFSigma);
					PSFy=kernel_IntGauss1D(jj, y, PSFSigma);
					
					kernel_DerivativeIntGauss1D(ii, x, PSFSigma, N, PSFy, &imdx, &imddx);
					kernel_DerivativeIntGauss1D(jj, y, PSFSigma, N, PSFx, &imdy, &imddy);

					//calculating second derivative og Intensity
					imddI =-1/pow(Nsigma, 2)/sz/sz;

					//denominator
					stepxtot  += imddx*cf-pow(imdx, 2)*df;
					stepytot  += imddy*cf-pow(imdy, 2)*df;
					stepItot  += imddI*cf-pow(imdI, 2)*df;

					//numerator
					dLLx+=imdx*cf;
					dLLy+=imdy*cf;
					dLLN+=imdI*cf;
				}

				x-=min(max(dLLx/stepxtot, -1e0f), 1e0f)/mm/2.0f;
				y-=min(max(dLLy/stepytot, -1e-0f), 1e-0f)/mm/2.0f;
				N-=min(max(dLLN/stepItot, -1e2f), 1e2f)/2.0f;
				N=max(N, 0.0001f);
				xarray[nn]=x;
				yarray[nn]=y;
				Narray[nn]=N;
			}
		}




		// calculating loglikelihood value

		Div=0.0f;
		xmin=1000;
		xmax=-1000;
		ymin=1000;
		ymax=-1000;

		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			// generating the model
			model=b;
			for (nn=0;nn<mm;nn++) {
				x=xarray[nn];
				y=yarray[nn];
				N=Narray[nn];

				if (x>xmax) xmax=x;
				if (x<xmin) xmin=x;
				if (y>ymax) ymax=y;
				if (y<ymin) ymin=y;


				PSFx=kernel_IntGauss1D(ii, x, PSFSigma);
				PSFy=kernel_IntGauss1D(jj, y, PSFSigma);
				model+=N*PSFx*PSFy;
			}

			data=s_data[sz*sz*tx+sz*jj+ii];

			/*if (model>0) {
			if (data>0)
			Div+=-2*(data*log(model)-model-data*log(data)+data);
			else
			Div+=-2*(-model);
			}*/
			if (data>0){
				Div+=-2*(data*log(model)-model-data*log(data)+data-0.5f*log(2*pi*data)-model*log(model)+model+model*log(model)-model+0.5f*log(2*pi*model));}

			else{
				Div+=-2*(-model-model*log(model)+model+model*log(model)-model+0.5f*log(2*pi*model));}
		}

		zval=sqrt(Div)-sqrt((float)(sz*sz-2*mm-1));
		pval=(zval<0)*(1-0.5f*exp(b1*zval+a1*pow(zval,2)))+(zval>0)*(0.5f*exp(b2*zval+a2*pow(zval,2)));
		if ((pval>maxpval) && (xmin>-round(1.5f*PSFSigma)) && (xmax <(sz-1+round(1.5f*PSFSigma))) && (ymin>-round(1.5f*PSFSigma)) && (ymax<(sz-1+round(1.5f*PSFSigma)))) {
			for (qq=0;qq<mm;qq++)

			{
				xarrayfit[qq]=xarray[qq];
				yarrayfit[qq]=yarray[qq];
				Narrayfit[qq]=Narray[qq];
			}
			maxpval=pval;
			minDiv=Div;
			bfit=b;
			nnfit=mm;

		}




		// judge by loglikelihood to either go through or not
		// Current setup only returns values that are good enough to pass
		// the threshold
	}

	//added 2011-10-15, intensity iteration

	Nsigma=0.0001;
	if ((xarrayfit[0]!=0)&&(maxpval>llThreshold*0.01)){
		b=bfit;
		for (qq=0;qq<nnfit;qq++)

		{
			xarray[qq]=xarrayfit[qq];
			yarray[qq]=yarrayfit[qq];
			Narray[qq]=Narrayfit[qq];

		}

		for (kk=0;kk<0;kk++){ //10 iteration for intensity estimation

			dLLb=0.0f;
			stepbgtot = 0;

			//generating the model and calc b update
			for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {

				model=b;
				for (nn=0;nn<nnfit;nn++) {
					x_current[nn]=xarray[nn];
					y_current[nn]=yarray[nn];
					N_current[nn]=Narray[nn];

					x=xarray[nn];
					y=yarray[nn];
					N=Narray[nn];

					PSFx=kernel_IntGauss1D(ii, x, PSFSigma);
					PSFy=kernel_IntGauss1D(jj, y, PSFSigma);
					model+=N*PSFx*PSFy;
				}
				data=s_data[sz*sz*tx+sz*jj+ii];
				cf=0.0f;
				df=0.0f;
				if (model>10e-3f) cf=data/model-1;
				if (model>10e-3f) df=data/pow(model, 2);
				cf=min(cf, 10e4f);
				df=min(df, 10e4f);
				stepbgtot += -df;
				dLLb+=cf;
			}
			b-=min(max(dLLb/stepbgtot, -1e0f), 1e0f)/nnfit/2;
			b=max(b, 0.001f);
			//b=0.01f;



			//This starts iterative routine for theta_i other than theta_bg which is calculated above.
			for (nn=0;nn<nnfit;nn++){


				dLLx=0.0f;
				dLLy=0.0f;
				dLLN=0.0f;
				stepxtot = 0.0f;
				stepytot = 0.0f;
				stepItot = 0.0f;

				for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
					//generate model using initial value or iteration value.
					data=s_data[sz*sz*tx+sz*jj+ii];
					model=b;
					for (nnM=0;nnM<nnfit;nnM++) {
						x=x_current[nnM];
						y=y_current[nnM];
						N=N_current[nnM];
						PSFx=kernel_IntGauss1D(ii, x, PSFSigma);
						PSFy=kernel_IntGauss1D(jj, y, PSFSigma);
						model+=N*PSFx*PSFy;
					}

					cf=0.0f;
					df=0.0f;
					if (model>10e-3f) cf=data/model-1;
					if (model>10e-3f) df=data/pow(model, 2);
					cf=min(cf, 10e4f);
					df=min(df, 10e4f);

					x=x_current[nn];
					y=y_current[nn];
					N=N_current[nn];

					PSFx=kernel_IntGauss1D(ii, x, PSFSigma);
					PSFy=kernel_IntGauss1D(jj, y, PSFSigma);
					
					kernel_DerivativeIntGauss1D(ii, x, PSFSigma, N, PSFy, &imdx, &imddx);
					kernel_DerivativeIntGauss1D(jj, y, PSFSigma, N, PSFx, &imdy, &imddy);

					// second Derivative of Intensity
					imddI =-1/pow(Nsigma, 2)/sz/sz;

					//denominator
					stepxtot  += imddx*cf-pow(imdx, 2)*df;
					stepytot  += imddy*cf-pow(imdy, 2)*df;
					stepItot  += imddI*cf-pow(imdI, 2)*df;

					//numerator
					dLLx+=imdx*cf;
					dLLy+=imdy*cf;
					dLLN+=imdI*cf;

				}

				x-=min(max(dLLx/stepxtot, -1e0f), 1e0f)/nnfit/2.0f;
				y-=min(max(dLLy/stepytot, -1e-0f), 1e-0f)/nnfit/2.0f;
				N-=min(max(dLLN/stepItot, -1e2f), 1e2f)/nnfit/2.0f;
				N=max(N, 0.0001f);
				xarray[nn]=x;
				yarray[nn]=y;
				Narray[nn]=N;
			}
		}


		//update LLR-pvalue after iteration of intensity

		// calculating loglikelihood value

		Div=0.0f;
		xmin=1000;
		xmax=-1000;
		ymin=1000;
		ymax=-1000;
		maxpval=-1.0f;
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			// generating the model
			model=b;
			for (nn=0;nn<nnfit;nn++) {
				x=xarray[nn];
				y=yarray[nn];
				N=Narray[nn];

				if (x>xmax) xmax=x;
				if (x<xmin) xmin=x;
				if (y>ymax) ymax=y;
				if (y<ymin) ymin=y;


				PSFx=kernel_IntGauss1D(ii, x, PSFSigma);
				PSFy=kernel_IntGauss1D(jj, y, PSFSigma);
				model+=N*PSFx*PSFy;
			}

			data=s_data[sz*sz*tx+sz*jj+ii];

			/*if (model>0) {
			if (data>0)
			Div+=-2*(data*log(model)-model-data*log(data)+data);
			else
			Div+=-2*(-model);
			}*/
			if (data>0){
				Div+=-2*(data*log(model)-model-data*log(data)+data-0.5f*log(2*pi*data)-model*log(model)+model+model*log(model)-model+0.5f*log(2*pi*model));}

			else{
				Div+=-2*(-model-model*log(model)+model+model*log(model)-model+0.5f*log(2*pi*model));}
		}

		zval=sqrt(Div)-sqrt((float)(sz*sz-2*nnfit-1));
		pval=(zval<0)*(1-0.5f*exp(b1*zval+a1*pow(zval,2)))+(zval>0)*(0.5f*exp(b2*zval+a2*pow(zval,2)));
		if ((xmin>-round(1.5f*PSFSigma)) && (xmax <(sz-1+round(1.5f*PSFSigma))) && (ymin>-round(1.5f*PSFSigma)) && (ymax<(sz-1+round(1.5f*PSFSigma)))) {
			for (qq=0;qq<nnfit;qq++)

			{
				xarrayfit[qq]=xarray[qq];
				yarrayfit[qq]=yarray[qq];
				Narrayfit[qq]=Narray[qq];
			}
			minDiv=Div;
			bfit=b;
			maxpval=pval;

		}
	}
	// output
	if (maxpval>llThreshold)
		goodness=1;

	int fitnum=num;

	for (nn=0;nn<fitnum;nn++)
		d_X[BlockSize*fitnum*bx+fitnum*tx+nn]=xarrayfit[nn]*goodness;

	for (nn=0;nn<fitnum;nn++)
		d_Y[BlockSize*fitnum*bx+fitnum*tx+nn]=yarrayfit[nn]*goodness;

	for (nn=0;nn<fitnum;nn++)
		d_N[BlockSize*fitnum*bx+fitnum*tx+nn]=Narrayfit[nn]*goodness;



	d_b[BlockSize*bx+tx]=bfit*goodness;
	d_Div[BlockSize*bx+tx]=minDiv;
	return;
}

//END OF KERNAL FUNCTION