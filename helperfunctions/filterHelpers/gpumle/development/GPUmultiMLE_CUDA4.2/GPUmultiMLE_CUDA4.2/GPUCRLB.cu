/*! \file GPUCRLB.cu
 *  \author Fang Huang
 *  \date October 10, 2010
 *  \brief This file contains the GPU code to perform Cramer-Rao Lower Bound
 *  calculations.
 */

#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "image_operation.h"

// Thread block size
#define BSZ 10
#define MATMEM 1240
#define MEM 60
#define IMSZ 11
#define IMSZBIG 21
#define NK 128 //number of blocks to run in each kernel
#define pi 3.141592f
#define NUM 5 //max number of multifit in previous step
#define min(a,b)            (((a) < (b)) ? (a) : (b))
//kernel_MLEFit<<<dimGrid, dimBlock>>>(ii, sz, BlockSize, fitnum, d_xarray, d_yarray, d_Narray, d_barray, d_fishermatrix, BlockSize);

void cudasafe( cudaError_t err, char* str, int lineNumber);
void CUDAERROR(const char *instr,int lineNumber);

__global__ void kernel_CRLB(int iiK,int blockx, int sz, int BlockSize, int fitnum,float Nave, 
							float Nsigma, float PSFSigma,float *d_xarray, float *d_yarray, 
							float *d_Narray,float *d_barray, float *d_CRLBarray , float *d_covariance); 


void GPUCRLB(int candicc,const int boxsz,float psfsigma,int fitnum, float *x, float *y, float *n, float *bg, float Nave,float Nsigmafit,float *CRLBarray, float *covariance)

{
/*!
 *  \brief Setup data for CUDA kernel_CRLB to compute CRLB values for the subregions.
 *  \param candicc
 *  \param boxsz
 *  \param psfsigma
 *  \param fitnum 
 *  \param x
 *	\param y
 *  \param n
 *  \param bg
 *  \param Nave
 *  \param Nsigmafit
 *  \param CRLBarray
 *  \param covariance
 */
	int blockx;
	int threadx;
	//int ii;
	float * d_xarray, *d_yarray, * d_Narray, * d_barray, * d_CRLBarray,*d_covariance;
	int framenum;
	//int Ndim;
	int  Nfits, totmatsize,memframenum;

	// declare text variable
	char text[1024];

	totmatsize=2*(2*fitnum+1)*(2*fitnum+1);
	int BlockSize =(int) floor((float)15000/4/totmatsize);
	BlockSize = max(4, BlockSize);
	BlockSize = min(BSZ, BlockSize);

	framenum=candicc;
	Nfits=framenum;

	memframenum=(int)ceil((float)framenum/BlockSize)+128;

	cudaMalloc((void**)&d_xarray, fitnum*memframenum*BlockSize*sizeof(float));
	cudaMemset(d_xarray, 0, fitnum*memframenum*BlockSize*sizeof(float));
	cudaMemcpy(d_xarray, x, fitnum*framenum*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_yarray, fitnum*memframenum*BlockSize*sizeof(float));
	cudaMemset(d_yarray, 0, fitnum*memframenum*BlockSize*sizeof(float));
	cudaMemcpy(d_yarray, y, fitnum*framenum*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_Narray, fitnum*memframenum*BlockSize*sizeof(float));
	cudaMemset(d_Narray, 0, fitnum*memframenum*BlockSize*sizeof(float));
	cudaMemcpy(d_Narray, n, fitnum*framenum*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_barray, memframenum*BlockSize*sizeof(float));
	cudaMemset(d_barray, 0, memframenum*BlockSize*sizeof(float));
	cudaMemcpy(d_barray, bg, framenum*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&d_CRLBarray, (2*fitnum+1)*memframenum*BlockSize*sizeof(float));
	cudaMemset(d_CRLBarray, 0, (2*fitnum+1)*memframenum*BlockSize*sizeof(float));

	cudaMalloc((void**)&d_covariance, fitnum*memframenum*BlockSize*sizeof(float));
	cudaMemset(d_covariance, 0, fitnum*memframenum*BlockSize*sizeof(float));

	//cudaMalloc((void**)&d_fishermatrix, (2*fitnum+1)*(2*fitnum+1)*memframenum*BlockSize*sizeof(float));
	//cudaMemset(d_fishermatrix, 0, (2*fitnum+1)*(2*fitnum+1)*memframenum*BlockSize*sizeof(float));

	int numK=(int)ceil((float)Nfits/BlockSize/NK);

	printf("Calculating CRLB...\n");
	for (int ii=0;ii<numK;ii++) {
		//blocksneed=((float) Nfits)/BlockSize;
		//setup kernel
		blockx = (int) min(ceil(((float)(((float)Nfits)/BlockSize)-ii*NK)), NK);
		blockx = max(blockx,1);
		threadx= BlockSize;


		dim3 dimBlock(threadx);
		dim3 dimGrid(blockx);

		



		kernel_CRLB<<<dimGrid, dimBlock>>>(ii,blockx,boxsz, BlockSize, fitnum,Nave,Nsigmafit, psfsigma,d_xarray, d_yarray, d_Narray, d_barray, d_CRLBarray,d_covariance);
		CUDAERROR(text,__LINE__);

		// too make the loop work, we have to operate on d_data
		// this trick works to over smart compiler!
		cudaMemcpy(d_xarray, x, 1*sizeof(float), cudaMemcpyHostToDevice);


		//printf("Calculating CRLB: %d/%d is DONE...\n", (ii+1), numK);
		mexEvalString("pause(0.001)");

	}

	cudaMemcpy(CRLBarray, d_CRLBarray, framenum*(2*fitnum+1)*sizeof(float), cudaMemcpyDeviceToHost);

	cudaMemcpy(covariance, d_covariance, framenum*fitnum*sizeof(float), cudaMemcpyDeviceToHost);

	cudaFree(d_xarray);
	cudaFree(d_yarray);
	cudaFree(d_Narray);
	cudaFree(d_barray);
	cudaFree(d_CRLBarray);
	cudaFree(d_covariance);

}


__global__ void kernel_CRLB(int iiK,int blockx, int sz, int BlockSize, int fitnum,float Nave, 
							float Nsigma, float PSFSigma,float *d_xarray, float *d_yarray, 
							float *d_Narray,float *d_barray, float *d_CRLBarray , float *d_covariance) 
{
/*!
 *  \brief Compute CRLB values for the specified subregions.
 *  \param iiK
 *  \param blockx
 *  \param fitnum
 *  \param sz 
 *  \param BlockSize
 *  \param fitnum
 *  \param Nave
 *  \param Nsigma
 *  \param PSFSigma
 *  \param d_xarray
 *  \param d_yarray
 *  \param d_Narray
 *  \param d_barray
 *  \param d_CRLBarray
 *  \param d_covariance
 */
	int tx = threadIdx.x; //matrix number index
	int bx = blockIdx.x;
	int kk, ll, ii, jj,nn,num;
	float x, y, N, b, tmp1;
	float PSFx, norm, PSFy, imdya, imdyb, imdxa, imdxb,imdx,imdy;
	//float imdI;
	float imdbg;
	int fitstate, matsz;
	int matsize=0,numco;
	float yy[17];
	float model,disx,disy;//
	__shared__ float s_xarray[MEM];
	__shared__ float s_yarray[MEM];
	__shared__ float s_Narray[MEM];
	__shared__ float s_barray[MEM];
	__shared__ float s_temp[MATMEM];
	__shared__ float s_fishermatrix[MATMEM];
	__shared__ float zerofisher[MATMEM];


	for (kk=0;kk<MATMEM;kk++){
		s_fishermatrix[kk]=0;
		s_temp[kk]=0;
	}
	norm=1.0f/2.0f/PSFSigma/PSFSigma;

	bx=bx+iiK*NK;
	//import datas from device to shared memory
	for (ii=0;ii<fitnum;ii++){
		s_xarray[fitnum*tx+ii]=d_xarray[fitnum*bx*BlockSize+fitnum*tx+ii];
		s_yarray[fitnum*tx+ii]=d_yarray[fitnum*bx*BlockSize+fitnum*tx+ii];
		s_Narray[fitnum*tx+ii]=d_Narray[fitnum*bx*BlockSize+fitnum*tx+ii];
	}

	s_barray[tx]=d_barray[bx*BlockSize+tx];

	//calculation starts.


	//    Div=0.0;Ixx=0.0;  Iyy=0.0;  III=0.0;  Ibgbg=0.0;Ixy=0.0;  IyI=0.0;  Iybg=0.0;
	//    IIbg=0.0; Ixbg=0.0; IxI=0.0;  imdy=0.0; imdx=0.0; imdI=0.0; imdbg=0.0;

	for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
		// generate model for pixel ii jj

		model=s_barray[tx];

		for (nn=0;nn<fitnum;nn++){
			fitstate=0;
			x=s_xarray[fitnum*tx+nn];
			y=s_yarray[fitnum*tx+nn];
			N=s_Narray[fitnum*tx+nn];
			if(N>10) fitstate=1;
			PSFx=1.0f/2.0f*(erf((ii-x+0.5f)*sqrt(norm))-erf((ii-x-0.5f)*sqrt(norm)))*fitstate;
			PSFy=1.0f/2.0f*(erf((jj-y+0.5f)*sqrt(norm))-erf((jj-y-0.5f)*sqrt(norm)))*fitstate;
			model+=N*PSFx*PSFy;
		}


		matsize=0;
		for (nn=0;nn<fitnum;nn++){
			fitstate=0;
			x=s_xarray[fitnum*tx+nn];
			y=s_yarray[fitnum*tx+nn];
			N=s_Narray[fitnum*tx+nn];
			if(N>10) fitstate=1;
			PSFx=1.0f/2.0f*(erf((ii-x+0.5f)*sqrt(norm))-erf((ii-x-0.5f)*sqrt(norm)))*fitstate;
			PSFy=1.0f/2.0f*(erf((jj-y+0.5f)*sqrt(norm))-erf((jj-y-0.5f)*sqrt(norm)))*fitstate;
			imdya = exp(-1.0f/2.0f*pow(((jj+0.5f-y)/PSFSigma), 2.0f))*fitstate;
			imdyb = exp(-1.0f/2.0f*pow(((jj-0.5f-y)/PSFSigma), 2.0f))*fitstate;
			imdxa = exp(-1.0f/2.0f*pow(((ii+0.5f-x)/PSFSigma), 2.0f))*fitstate;
			imdxb = exp(-1.0f/2.0f*pow(((ii-0.5f-x)/PSFSigma), 2.0f))*fitstate;

			//calculating first derivatives

			imdx    = -N/sqrt(2.0f*pi)/PSFSigma*(imdxa-imdxb)*PSFy;
			imdy    = -N/sqrt(2.0f*pi)/PSFSigma*(imdya-imdyb)*PSFx;
			//imdI    = PSFx*PSFy-1/pow(Nsigma,2)*(N-Nave)/sz/sz*fitstate;

			//store in temp memory (which's allocate for LUDC method)
			s_temp[(2*fitnum+1)*tx+2*nn]=imdx;
			s_temp[(2*fitnum+1)*tx+2*nn+1]=imdy;
			//s_temp[(3*fitnum+1)*tx+3*nn+2]=imdI;
			matsize+=fitstate;
		}
		imdbg = 1.0f;
		// put derivitive of bg right after fitted parameters
		if (matsize>0)
			s_temp[(2*fitnum+1)*tx+2*matsize]=imdbg;


		for (kk=0;kk<(2*matsize+1);kk++)
			for (ll=0;ll<(2*matsize+1);ll++)
				s_fishermatrix[tx*(2*fitnum+1)*(2*fitnum+1)+kk*(2*matsize+1)+ll]+=s_temp[(2*fitnum+1)*tx+kk]*s_temp[(2*fitnum+1)*tx+ll]/model;

	}



	for (kk=0;kk<MATMEM;kk++){

		s_temp[kk]=0;
	}
	/*-----------------------------------------------------------------------------------------------------------------------
	-----------------------------------------------------------------------------------------------------------------------*/
	//cope fisher information matrix into device memory for output.

	//for (kk=0;kk<(2*matsize+1);kk++)
	//	for (jj=0;jj<(2*matsize+1);jj++)
	//		d_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*bx*BlockSize+(2*fitnum+1)*(2*fitnum+1)*tx+kk*(2*fitnum+1)+jj]=s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*tx+kk*(2*matsize+1)+jj];
	//modification of fisher information matrix
	matsz=matsize*2+1;
	if(matsize>0)
	{
		//copy fisher matrix
		for (ii=0;ii<matsz;ii++)
			for(kk=0;kk<matsz;kk++)
				zerofisher[(2*fitnum+1)*(2*fitnum+1)*tx+kk*(2*matsize+1)+ii]=s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*tx+kk*(2*matsize+1)+ii];

		//set covariance to zero
		for(kk=0;kk<matsize;kk++){
			for(jj=0;jj<matsize;jj++)
			{
				if(kk!=jj)
				{
					zerofisher[(2*fitnum+1)*(2*fitnum+1)*tx+(2*kk+1)*(2*matsize+1)+(2*jj+1)]=0;
					zerofisher[(2*fitnum+1)*(2*fitnum+1)*tx+(2*kk)*(2*matsize+1)+(2*jj)]=0;
				}
			}
		}


		//inverse 0 fisher

		tmp1=0;
		for (jj = 0; jj < matsz; jj++) {
			//calculate upper matrix
			for (ii=0;ii<=jj;ii++)   {
				//deal with ii-1 in the sum, set sum(kk=0->ii-1) when ii=0 to zero
				if (ii>0) {
					for (kk=0;kk<=ii-1;kk++) {
						tmp1=tmp1+zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+kk*matsz]*zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+kk+jj*matsz];
					}
					zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz]=zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz]-tmp1;
					tmp1=0;
				}
				else {

					zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz]=zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz];
					zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz]=zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz];

				}
			}


			//calculate lowwer matrix
			for (ii=jj+1;ii<matsz;ii++) {
				if (jj>0) {
					for (kk=0;kk<=jj-1;kk++) {
						tmp1=tmp1+zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+kk*matsz]*zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+kk+jj*matsz];
					}

					zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz]=(1/zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+jj+jj*matsz])*(zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz]-tmp1);

					tmp1=0;
				}
				else {
					zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz]=(1/zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+jj+jj*matsz])*zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz];
					zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz]=zerofisher[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz];
				}
			}

		}





		tmp1=0;

		for (num=0;num<matsz;num++) {
			// calculate yy
			if (num==0)
				yy[0]=1;

			else
				yy[0]=0;

			for (ii=1;ii<matsz;ii++) {

				if (ii==num)
					b=1;
				else
					b=0;

				for (jj=0;jj<=ii-1;jj++) {
					tmp1=zerofisher[threadIdx.x*(2*fitnum+1)*(2*fitnum+1)+ii+jj*matsz]*yy[jj]+tmp1;
				}
				yy[ii]=b-tmp1;
				tmp1=0;

			}


			// calculate xx
			s_temp[threadIdx.x*(2*fitnum+1)*(2*fitnum+1)+matsz-1+num*matsz]=yy[matsz-1]/zerofisher[threadIdx.x*(2*fitnum+1)*(2*fitnum+1)+(matsz-1)+(matsz-1)*matsz];

			for (ii=matsz-2;ii>=0;ii--) {
				for (jj=ii+1;jj<matsz;jj++) {
					tmp1=zerofisher[threadIdx.x*(2*fitnum+1)*(2*fitnum+1)+ii+jj*matsz]*s_temp[threadIdx.x*(2*fitnum+1)*(2*fitnum+1)+jj+num*matsz]+tmp1;
				}
				s_temp[threadIdx.x*(2*fitnum+1)*(2*fitnum+1)+ii+num*matsz]=(1/zerofisher[threadIdx.x*(2*fitnum+1)*(2*fitnum+1)+ii+ii*matsz])*(yy[ii]-tmp1);

				tmp1=0;
			}
		}



		//Calculate new fisher
		for(kk=0;kk<matsize;kk++)
		{
			s_Narray[tx*fitnum+kk]=sqrt(abs(s_temp[tx*(2*fitnum+1)*(2*fitnum+1)+2*kk*matsz+2*kk]));        //crlbx
			s_barray[tx*fitnum+kk]=sqrt(abs(s_temp[tx*(2*fitnum+1)*(2*fitnum+1)+(2*kk+1)*matsz+(2*kk+1)])); //crlby
		}


		for(kk=0;kk<matsize;kk++)
		{
			for(jj=0;jj<matsize;jj++)
			{
				if(kk!=jj)
				{
					disx=pow((abs(s_xarray[tx*fitnum+kk]-s_xarray[tx*fitnum+jj])),2)/s_Narray[tx*fitnum+kk]/s_Narray[tx*fitnum+jj];
					disy=pow((abs(s_yarray[tx*fitnum+kk]-s_yarray[tx*fitnum+jj])),2)/s_barray[tx*fitnum+kk]/s_barray[tx*fitnum+jj];
					s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*tx+(2*kk+1)*(2*matsize+1)+(2*jj+1)]=disy/(disy+1)*s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*tx+(2*kk+1)*(2*matsize+1)+(2*jj+1)];
					s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*tx+(2*kk)*(2*matsize+1)+(2*jj)]=disx/(disx+1)*s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*tx+(2*kk)*(2*matsize+1)+(2*jj)];
				}
			}

		}

		//inverse new fisher
		// matsize will be a indicator saying how big the matrix is for each thread. matsize*2+1 will be the scale of the matrix.
		// start calculating matrix inverse
		/*-----------------------------------------------------------------------------------------------------------------------
		-----------------------------------------------------------------------------------------------------------------------*/
		for (kk=0;kk<MATMEM;kk++){

			s_temp[kk]=0;
		}

		tmp1=0;
		for (jj = 0; jj < matsz; jj++) {
			//calculate upper matrix
			for (ii=0;ii<=jj;ii++)   {
				//deal with ii-1 in the sum, set sum(kk=0->ii-1) when ii=0 to zero
				if (ii>0) {
					for (kk=0;kk<=ii-1;kk++) {
						tmp1=tmp1+s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+kk*matsz]*s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+kk+jj*matsz];
					}
					s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz]=s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz]-tmp1;
					tmp1=0;
				}
				else {

					s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz]=s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz];
					s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz]=s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz];

				}
			}


			//calculate lowwer matrix
			for (ii=jj+1;ii<matsz;ii++) {
				if (jj>0) {
					for (kk=0;kk<=jj-1;kk++) {
						tmp1=tmp1+s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+kk*matsz]*s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+kk+jj*matsz];
					}

					s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz]=(1/s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+jj+jj*matsz])*(s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz]-tmp1);

					tmp1=0;
				}
				else {
					s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz]=(1/s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+jj+jj*matsz])*s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz];
					s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz]=s_fishermatrix[(2*fitnum+1)*(2*fitnum+1)*threadIdx.x+ii+jj*matsz];
				}
			}

		}





		tmp1=0;

		for (num=0;num<matsz;num++) {
			// calculate yy
			if (num==0)
				yy[0]=1;

			else
				yy[0]=0;

			for (ii=1;ii<matsz;ii++) {

				if (ii==num)
					b=1;
				else
					b=0;

				for (jj=0;jj<=ii-1;jj++) {
					tmp1=s_fishermatrix[threadIdx.x*(2*fitnum+1)*(2*fitnum+1)+ii+jj*matsz]*yy[jj]+tmp1;
				}
				yy[ii]=b-tmp1;
				tmp1=0;

			}


			// calculate xx
			s_temp[threadIdx.x*(2*fitnum+1)*(2*fitnum+1)+matsz-1+num*matsz]=yy[matsz-1]/s_fishermatrix[threadIdx.x*(2*fitnum+1)*(2*fitnum+1)+(matsz-1)+(matsz-1)*matsz];

			for (ii=matsz-2;ii>=0;ii--) {
				for (jj=ii+1;jj<matsz;jj++) {
					tmp1=s_fishermatrix[threadIdx.x*(2*fitnum+1)*(2*fitnum+1)+ii+jj*matsz]*s_temp[threadIdx.x*(2*fitnum+1)*(2*fitnum+1)+jj+num*matsz]+tmp1;
				}
				s_temp[threadIdx.x*(2*fitnum+1)*(2*fitnum+1)+ii+num*matsz]=(1/s_fishermatrix[threadIdx.x*(2*fitnum+1)*(2*fitnum+1)+ii+ii*matsz])*(yy[ii]-tmp1);

				tmp1=0;
			}
		}

	}






	// finished    
	// copy back to device memory
	numco=0;
	for (kk=0;kk<matsz;kk++)
	{
		d_CRLBarray[(2*fitnum+1)*bx*BlockSize+tx*(2*fitnum+1)+kk]=sqrt(abs(s_temp[tx*(2*fitnum+1)*(2*fitnum+1)+kk+kk*matsz]));
		if ((kk+1)%2==0)
		{
			d_covariance[fitnum*bx*BlockSize+tx*fitnum+numco]=s_temp[tx*(2*fitnum+1)*(2*fitnum+1)+kk+(kk-1)*matsz];
			numco++;
		}
	}


	return;



}


