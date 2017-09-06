/*! \file fitmultimain.cu
 *  \author Fang Huang
 *  \date October 10, 2010
 *  \brief This file contains the mexFunction() call and relevant matlab
 *  interface code.
 */
/*
*float  	Floating point number.  	4bytes  	+/- 3.4e +/- 38 (~7 digits)
*/
#include <windows.h>
#pragma comment(lib, "kernel32.lib")

#include <math.h>
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include "image_operation.h"
#include "filter.h"
#include "definitions.h"

void maxfilter(const float * data , float * result, const int fsz , const int xsz, const int ysz, const int tsz); 

void cudasafe( cudaError_t err, char* str, int lineNumber);
void CUDAERROR(const char *instr,int lineNumber);

void CUDAERROR(const char *instr,int lineNumber) {
	cudaError_t errornum;
	const char *str;
	if (errornum=cudaGetLastError()) {
		//reset all cuda devices
		int deviceCount = 0;
		int ii = 0;
		cudasafe(cudaGetDeviceCount(&deviceCount),"cudaGetDeviceCount",__LINE__ ); //query number of GPUs
		for (ii = 0; ii< deviceCount;ii++) {
			cudaSetDevice(ii);
			cudaDeviceReset();
		}
		str=cudaGetErrorString(errornum);
		//mexPrintf("gpuGaussNDmle(line %i): %s in %s\n",lineNumber, str, instr);
		cudaDeviceReset();
		mexErrMsgIdAndTxt("gpuGaussNDmle:cudaFail","gpuGaussNDmle(line %i): %s in %s\n",lineNumber, str, instr);
		exit(1); // might not stop matlab
	}
}

void cudasafe( cudaError_t err, char* str, int lineNumber)
{
	if (err != cudaSuccess)
	{
		//reset all cuda devices
		int deviceCount = 0;
		int ii = 0;
		cudasafe(cudaGetDeviceCount(&deviceCount),"cudaGetDeviceCount",__LINE__ ); //query number of GPUs
		for (ii = 0; ii< deviceCount;ii++) {
			cudaSetDevice(ii);
			cudaDeviceReset();
		}
		mexErrMsgIdAndTxt("gpuGaussNDmle:cudaFail","%s failed with error code %i at line %d\n",str,err, lineNumber);
		exit(1); // might not stop matlab
	}
}

bool mxIsScalar(const mxArray *array_ptr)
{
/*!
 *  \brief Check that the passed in mxArray is infact a single number.
 *  \param array_ptr pointer to the array to check.
 *	\return bool
 */
	const mwSize *size;
	const size_t dims = mxGetNumberOfDimensions(array_ptr);
	size = mxGetDimensions(array_ptr);
	if (dims != 2) return false;
	if (size[0] != 1) return false;
	if (size[1] != 1) return false;
	return true;
}

void printTiming(char *var, LARGE_INTEGER start) {
/*!
 *  \brief This returns the runtime information to Matlab using the high precision timer.
 *  \param var Name of the variable to assign the timing info to.
 *  \param start The starting time for this event.
 */
	LARGE_INTEGER freq, stop;
	char timing[100];
	QueryPerformanceFrequency( &freq );
	QueryPerformanceCounter(&stop);
	sprintf_s(timing,(size_t) 100,"%s=%f", var, (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart);
	mexEvalString(timing);
}

void mexFunction(int nlhs, mxArray *plhs[],	int	nrhs, const	mxArray	*prhs[])

{
/*!
 *  \brief Entry point in the code for Matlab.  Equivalent to main().
 *  \param nlhs number of left hand mxArrays to return
 *  \param plhs array of pointers to the output mxArrays
 *  \param nrhs number of input mxArrays
 *  \param prhs array of pointers to the input mxArrays.
 */
	int ii=0,jj=0,kk=0,tt=0,candicc=0,candi=0;  //!< various counters unsed throughout the code
	const mwSize *pSizeA=0;
	mwSize outsize[2],outsize2[2];
	const float *indata=0;
	int boxsz=0;
	float psfsigma=0;
	int fitnum=0;
	float Nave=0;
	float *unif1=0,*unif2=0,*unif=0;  //!< uniformed filtered copies of the subregion
	float *maxf=0;					  //!< maximum filtered copy of the subregion
	float *candx=0,*candy=0,*candz=0,*left=0,*top=0,*in_sub=0;
	//float *out;

	float threshold;
	size_t Adims,data_bytes;
	float *xx,*yy,*nn,*x,*y,*n,*bb,*div,testx,testy,testn;
	float *CRLBarray,*covariance;
	float *SuperResim;
	float resolution;
	float pixelsz;
	float zoom;
	float Nsigmafit=0;
	//fun with timing
	LARGE_INTEGER t1, t2, t3, t4;

	printf("This code developed development by Lidke Lab at UNM. \nThe author(Fang Huang, Keith Lidke, and Lidke Lab's other group member)\nreserve all rights of using this code.\nReferece: Simultaneous multiple-emitter fitting for single molecule super-resolution imaging, \nFang Huang, Samantha L. Schwartz, Jason M. Byars, and Keith A. Lidke,\nBiomedical Optics Express, Vol. 2, Issue 5, pp. 1377-1393 (2011)\n" );

	printf("Start Fitting...\n" );

	//mexEvalString("pause(0.1)");

	if (nrhs != 9) {//																			0											1					2						3				4				5			       6	7	8	
		mexErrMsgIdAndTxt("MATLAB:WrongNumberOfInputs", "This function needs 10 inputs: 3D-image(x,y,frame,matlab matrix,single),PSFsigma(in pixels),initial estimated intensity,fit_type(1-5),pvaluethreshold (0.05),resolution(nm),pixelsz(nm),zm,boxsz");
	}


	Adims=mxGetNumberOfDimensions(prhs[0]);

	if (mxGetClassID(prhs[0])!=mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("Multifit:WrongPrecision","3D-Data must be comprised of single floats! Please use single(data) to convert into correct input data type\n");

	if (Adims != 3)
		mexErrMsgIdAndTxt("Multifit:WrongDimensions", "Only 3D data set can be fit by multi fluorophore fitting. X by Y by Frame");

	pSizeA = mxGetDimensions(prhs[0]);
	const int data_xsz=(int) pSizeA[0];
	const int data_ysz=(int) pSizeA[1];
	const int data_frame=(int) pSizeA[2];
	if (data_xsz <16 || data_ysz <16)
		mexErrMsgTxt("3D-Data must be at least 16 by 16 by frame for this fitting.");

	if (!mxIsNumeric(prhs[1]))
		mexErrMsgIdAndTxt("Multifit:NotNumeric", "psfsigma must be a numeric value\n");
	if (!mxIsNumeric(prhs[2]))
		mexErrMsgIdAndTxt("Multifit:NotNumeric", "estimated must be a numeric value\n");
	if (!mxIsNumeric(prhs[3]))
		mexErrMsgIdAndTxt("Multifit:NotNumeric", "fit_type must be a numeric value\n");
	if (!mxIsNumeric(prhs[4]))
		mexErrMsgIdAndTxt("Multifit:NotNumeric", "threshold must be a numeric value\n");
	if (!mxIsNumeric(prhs[5]))
		mexErrMsgIdAndTxt("Multifit:NotNumeric", "resolution must be a numeric value\n");
	if (!mxIsNumeric(prhs[6]))
		mexErrMsgIdAndTxt("Multifit:NotNumeric", "pixel must be a numeric value\n");
	if (!mxIsNumeric(prhs[7]))
		mexErrMsgIdAndTxt("Multifit:NotNumeric", "zoom must be a numeric value\n");

	if (!mxIsScalar(prhs[1]))
		mexErrMsgIdAndTxt("Multifit:WrongDimensions","psfsigma must be scalar\n");
	if (!mxIsScalar(prhs[2]))
		mexErrMsgIdAndTxt("Multifit:WrongDimensions","estimated intensity must be scalar\n");
	if (!mxIsScalar(prhs[3]))
		mexErrMsgIdAndTxt("Multifit:WrongDimensions","fit_type must be scalar\n");
	if (!mxIsScalar(prhs[4]))
		mexErrMsgIdAndTxt("Multifit:WrongDimensions","threshold must be scalar\n");
	if (!mxIsScalar(prhs[5]))
		mexErrMsgIdAndTxt("Multifit:WrongDimensions","resolution must be scalar\n");
	if (!mxIsScalar(prhs[6]))
		mexErrMsgIdAndTxt("Multifit:WrongDimensions","pixel size must be scalar\n");
	if (!mxIsScalar(prhs[7]))
		mexErrMsgIdAndTxt("Multifit:WrongDimensions","zoom must be scalar\n");

	indata=(float *)mxGetData(prhs[0]);
	psfsigma=(float)mxGetScalar(prhs[1]);
	Nave=(float)mxGetScalar(prhs[2]);
	fitnum=(int)mxGetScalar(prhs[3]);
	threshold=(float)mxGetScalar(prhs[4]);
	resolution=(float)mxGetScalar(prhs[5]);
	pixelsz=(float)mxGetScalar(prhs[6]);
	zoom=(float)mxGetScalar(prhs[7]);
	boxsz=(int)mxGetScalar(prhs[8]);
	//Nsigmafit=(float)mxGetScalar(prhs[9]);
	Nsigmafit=0.00001; // fixed intensity for position initial estimates, this restriction is lifted on 2nd round iteration of itensity,position fitting.
	if (psfsigma<=0.2f) {
		mexErrMsgIdAndTxt("Multifit:InputOutofRange","PSFsigma must be greater than 0.2.\nThe camera pixel size is too large." );
	}
	if (psfsigma>4.0f) {
		mexErrMsgIdAndTxt("Multifit:InputOutofRange","PSFsigma should be less than 4.0.\nPlease check your sigma value." );
	}

	if ((fitnum<1)||(fitnum>8))
		mexErrMsgTxt("Fitting model number out of range. Fitnum has to be within the range of 1-8.");

	//if (threshold<=0 || threshold > 1000)
	//	mexErrMsgIdAndTxt("Multifit:InputOutofRange","NLLR threshold is out of range. Please use any number between 0.0001-1000");

	if (resolution<=0) {
		resolution=20;
		mexErrMsgIdAndTxt("Multifit:InputOutofRange","resolution must be greater than 0. Default value is 20(nm).\n" );

	}

	if (pixelsz<=0) {
		pixelsz=106;
		mexErrMsgIdAndTxt("Multifit:InputOutofRange","Pixel Size must be greater than 0. Default value is 106(nm).\n" );

	}

	//What is smallest possible value?
	if (zoom<=0) {
		zoom=20;
		mexErrMsgIdAndTxt("Multifit:InputOutofRange","Zoom must be greater than 0. Default value is 30.\n" );

	}

	if (zoom - floor(zoom) > 0.000001f) {
		zoom = ceil(zoom);
		mexWarnMsgIdAndTxt("Multifit:InputOutofRange","Zoom should be an integer. Rounding up to next whole number value.\n");
	}

	data_bytes=data_xsz*data_ysz*data_frame*sizeof(float);

	indata=(float *)mxGetData(prhs[0]);

	data_bytes=((float) pSizeA[0]*pSizeA[1]*pSizeA[2])*sizeof(float);

	//uniform and max filter is used to find candidates in order to cut box from.
	//uniform max filter to find subregions
	const int unifsz=(int) round(psfsigma*2+1);
	const int maxfsz=(int) (boxsz-1);

	unif1=(float*) malloc (data_bytes); //need to be freed.
	memset(unif1,0,data_bytes);
	unif2=(float*) malloc (data_bytes); //need to be freed.
	memset(unif2,0,data_bytes);
	unif=(float*) malloc (data_bytes); //need to be freed.
	memset(unif,0,data_bytes);
	maxf=(float*) malloc (data_bytes); //need to be freed.
	memset(maxf,0,data_bytes);


	QueryPerformanceCounter(&t1);
	unifilter(indata , unif1 , unifsz,data_xsz,data_ysz,data_frame);
	unifilter(indata , unif2 , 2*unifsz,data_xsz,data_ysz,data_frame);
	
	
	//int off = 0;
	/*for(tt=0;tt<data_frame;tt++) {
		for(jj=0;jj<data_ysz;jj++) {
			for(ii=0;ii<data_xsz;ii++) {
				off = tt*data_xsz*data_ysz+data_xsz*jj+ii;
				unif[tt*data_xsz*data_ysz+data_xsz*jj+ii]=unif1[off]-unif2[off];
			}
		}
	}*/
	for (ii=0; ii<data_frame*data_ysz*data_xsz; ii++)
		unif[ii]=unif1[ii]-unif2[ii];

	free(unif1);
	free(unif2);
	maxfilter(unif , maxf , maxfsz,data_xsz,data_ysz,data_frame);

	float minI=Nave/2/pi/2/psfsigma/psfsigma/3;
	//minI=25;
	candicc=0;
	//find number of candidates for malloc
	/*for(tt=0;tt<data_frame;tt++)
		for(ii=0;ii<data_xsz;ii++)
			for(jj=0;jj<data_ysz;jj++)
	{
		kk=tt*data_xsz*data_ysz+data_xsz*jj+ii;
	
		if (((0.999*maxf[kk])<=unif[kk]) && (unif[kk]>=minI))
			candicc++;
	
	}*/
	for (ii=0; ii<data_frame*data_ysz*data_xsz; ii++)
		if (((0.999*maxf[ii])<=unif[ii]) && (unif[ii]>=minI))
			candicc++;
	
	//alloc memory
	candx=(float*) malloc (candicc*sizeof(float));//cand x, need to be freed.
	memset(candx,0,candicc*sizeof(float));
	candy=(float*) malloc (candicc*sizeof(float));//cand y, need to be freed.
	memset(candy,0,candicc*sizeof(float));
	candz=(float*) malloc (candicc*sizeof(float));//cand z, need to be freed.
	memset(candz,0,candicc*sizeof(float));

	candi=0;
	//store candidates
	for(tt=0;tt<data_frame;tt++)
		for(jj=0;jj<data_ysz;jj++)
			for(ii=0;ii<data_xsz;ii++)
	{
		kk=tt*data_xsz*data_ysz+data_xsz*jj+ii;

		if (((0.999f*maxf[kk])<=unif[kk]) && (unif[kk]>=minI))
		{
			candx[candi]=(float) ii;//try
			candy[candi]=(float) jj;
			candz[candi]=(float) tt;
			candi++;
		}

	}
	free(unif);
	free(maxf);


	// Cut subregions around found candidates
	//const int boxsz=5*psfsigma+1;
	left=(float*) malloc (candicc*sizeof(float));
	memset(left,0,candicc*sizeof(float));// left for cut sub regions, need to be freed
	top=(float*) malloc (candicc*sizeof(float)); // top for cut sub regions, need to be freed
	memset(top,0,candicc*sizeof(float));
	in_sub=(float*) malloc (boxsz*boxsz*candicc*sizeof(float));
	memset(in_sub,0,boxsz*boxsz*candicc*sizeof(float));

	//test
	/* outsize2[0]=boxsz;
	outsize2[1]=boxsz;
	outsize2[2]=candicc;
	plhs[1]= mxCreateNumericArray(Adims ,outsize2, mxSINGLE_CLASS, mxREAL);
	out=(float *)mxGetData(plhs[1]); */



	MakeSubregions(candicc, indata, data_xsz, data_ysz, data_frame, candx,candy, candz, (float) boxsz, in_sub,left,top);
	free(candx);
	free(candy);


	printTiming("t1",t1);



	//Fitting started
	x=(float*) malloc (fitnum*candicc*sizeof(float));//x, need to be freed.
	memset(x,0,fitnum*candicc*sizeof(float));
	y=(float*) malloc (fitnum*candicc*sizeof(float));//y, need to be freed.
	memset(y,0,fitnum*candicc*sizeof(float));
	n=(float*) malloc (fitnum*candicc*sizeof(float));//y, need to be freed.
	memset(n,0,fitnum*candicc*sizeof(float));
	bb=(float*) malloc (candicc*sizeof(float));//bg, need to be freed.
	memset(bb,0,candicc*sizeof(float));
	div=(float*) malloc (candicc*sizeof(float));//bg, need to be freed.
	memset(div,0,candicc*sizeof(float));
	//const float Nsigmafit=0.0001;
	const int iterations=50;

	/*float *out3;
	mwSize outsize3[2];
	outsize3[0]=candicc;
	outsize3[1]=1;
	//outsize3[2]=candicc;
	plhs[3]= mxCreateNumericArray(2 ,outsize3, mxSINGLE_CLASS, mxREAL);
	out3=(float *)mxGetData(plhs[3]);

	*/
	QueryPerformanceCounter(&t2);
	GPUmultifit(candicc,in_sub,boxsz,psfsigma,iterations,fitnum,Nave,Nsigmafit,threshold,x,y,n,bb,div);


	printTiming("t2",t2);
	QueryPerformanceCounter(&t3);

	mexEvalString("s3=tic");

	float *bbb,testb;
	// Boundary limitation
	xx=(float*) malloc (fitnum*candicc*sizeof(float));//x, need to be freed.
	memset(xx,0,fitnum*candicc*sizeof(float));
	yy=(float*) malloc (fitnum*candicc*sizeof(float));//y, need to be freed.
	memset(yy,0,fitnum*candicc*sizeof(float));
	nn=(float*) malloc (fitnum*candicc*sizeof(float));//y, need to be freed.
	memset(nn,0,fitnum*candicc*sizeof(float));
	bbb=(float*) malloc (fitnum*candicc*sizeof(float));//y, need to be freed.
	memset(bbb,0,fitnum*candicc*sizeof(float));




	/*float *out3;
	mwSize outsize3[2];
	outsize3[0]=fitnum;
	outsize3[1]=candicc;
	//outsize3[2]=candicc;
	plhs[3]= mxCreateNumericArray(2 ,outsize3, mxSINGLE_CLASS, mxREAL);
	out3=(float *)mxGetData(plhs[3]);
	*/

	int *fittot;
	fittot=(int*) malloc(candicc*sizeof(int));
	memset(fittot,0,candicc*sizeof(int));

	for(ii=0;ii<candicc;ii++)
	{
		kk=0;
		for(jj=0;jj<fitnum;jj++)
		{
			if ( (abs(x[fitnum*ii+jj])>0.001) && (abs(x[fitnum*ii+jj])<100) )
				fittot[ii]++;
			//out3[ii*fitnum+jj]=x[ii*fitnum+jj];//testx
			testx=x[ii*fitnum+jj];
			testy=y[ii*fitnum+jj];
			testn=n[ii*fitnum+jj];
			testb=bb[ii];
			if (true)
			{
				xx[ii*fitnum+kk]=testx;//testx
				yy[ii*fitnum+kk]=testy;
				nn[ii*fitnum+kk]=testn;
				bbb[ii*fitnum+kk]=testb;
				kk++;
			}
		}
	}


	free(x);
	free(y);
	free(n);
	free(bb);

	//Calculate CRLB

	CRLBarray=(float*) malloc ((fitnum*2+1)*candicc*sizeof(float));//x, need to be freed.
	memset(CRLBarray,0,(fitnum*2+1)*candicc*sizeof(float));
	covariance=(float*) malloc (fitnum*candicc*sizeof(float));//y, need to be freed.
	memset(covariance,0,fitnum*candicc*sizeof(float));


	//test




	GPUCRLB(candicc, boxsz, psfsigma,fitnum, xx, yy, nn, bbb, Nave, Nsigmafit,CRLBarray,covariance);


	int *fittmp;
	fittmp=(int*) malloc(candicc*sizeof(int));
	memset(fittmp,0,candicc*sizeof(int));

	float LAlim=resolution/pixelsz/1.7f;
	// filtering and assembling according to CRLB
	int recontot=0;
	for(ii=0;ii<candicc;ii++){
		for(jj=0;jj<fitnum;jj++){
			if ( (abs(xx[fitnum*ii+jj])>0.001) && (abs(xx[fitnum*ii+jj])<100) ) {
				fittmp[ii]++;
			}
			if ((xx[fitnum*ii+jj]>0.01f)&&(xx[fitnum*ii+jj]<(boxsz-1.01f))&&
				(yy[fitnum*ii+jj]>0.01f)&&(yy[fitnum*ii+jj]<(boxsz-1.01f))&&
				(nn[fitnum*ii+jj]>0)&&(CRLBarray[(fitnum*2+1)*ii+jj]<LAlim)&&
				(CRLBarray[(fitnum*2+1)*ii+jj]>0.0001f)&&
				(CRLBarray[(fitnum*2+1)*ii+jj+1]<LAlim)&&
				(CRLBarray[(fitnum*2+1)*ii+jj+1]>0.0001f))
				recontot++;
		}
	}

	
	float *xbf=(float*) malloc(recontot*sizeof(float));
	memset(xbf,0,recontot*sizeof(float));
	float *bbf=(float*) malloc(recontot*sizeof(float));
	memset(bbf,0,recontot*sizeof(float));
	float *ybf=(float*) malloc(recontot*sizeof(float));
	memset(ybf,0,recontot*sizeof(float));
	float *nbf=(float*) malloc(recontot*sizeof(float));
	memset(nbf,0,recontot*sizeof(float));
	float *tbf=(float*) malloc(recontot*sizeof(float));
	memset(tbf,0,recontot*sizeof(float));
	float *uncerxbf=(float*) malloc(recontot*sizeof(float));
	memset(uncerxbf,0,recontot*sizeof(float));
	float *uncerybf=(float*) malloc(recontot*sizeof(float));
	memset(uncerybf,0,recontot*sizeof(float));
	float *uncerbgbf=(float*) malloc(recontot*sizeof(float));
	memset(uncerbgbf,0,recontot*sizeof(float));
	float *covbf=(float*) malloc(recontot*sizeof(float));
	memset(covbf,0,recontot*sizeof(float));
	float *NLLRbf=(float*) malloc(recontot*sizeof(float));
	memset(NLLRbf,0,recontot*sizeof(float));
	float *fitNbf=(float*) malloc(recontot*sizeof(float));
	memset(fitNbf,0,recontot*sizeof(float));
	float *topbf=(float*) malloc(recontot*sizeof(float));
	memset(topbf,0,recontot*sizeof(float));
	float *leftbf=(float*) malloc(recontot*sizeof(float));
	memset(leftbf,0,recontot*sizeof(float));
	/*	outsize2[0]=2*fitnum+1;
	outsize2[1]=candicc;
	plhs[1]= mxCreateNumericArray(2 ,outsize2, mxSINGLE_CLASS, mxREAL);
	out=(float *)mxGetData(plhs[1]); */

	//adding x,y,t,crlbx,crlby,crlb-bg,NLLR.



	//test
	//float *cov;
	//plhs[4]= mxCreateNumericArray(2 ,outsize2, mxSINGLE_CLASS, mxREAL);
	//cov=(float *)mxGetData(plhs[4]);
	//test-end

	// before filtering coords collections.
	int reconii=0;
	for(ii=0;ii<candicc;ii++){
		for(jj=0;jj<fitnum;jj++){
			if ((xx[fitnum*ii+jj]>0.01f)&&(xx[fitnum*ii+jj]<(boxsz-1.01f))&&
				(yy[fitnum*ii+jj]>0.01f)&&(yy[fitnum*ii+jj]<(boxsz-1.01f))&&
				(nn[fitnum*ii+jj]>0)&&(CRLBarray[(fitnum*2+1)*ii+jj]<LAlim)&&
				(CRLBarray[(fitnum*2+1)*ii+jj]>0.0001f)&&
				(CRLBarray[(fitnum*2+1)*ii+jj+1]<LAlim)&&
				(CRLBarray[(fitnum*2+1)*ii+jj+1]>0.0001f)
				&&(covariance[fitnum*ii+jj]/CRLBarray[(fitnum*2+1)*ii+jj]/CRLBarray[(fitnum*2+1)*ii+jj+1]<1)
				&&(covariance[fitnum*ii+jj]/CRLBarray[(fitnum*2+1)*ii+jj]/CRLBarray[(fitnum*2+1)*ii+jj+1]>-1))
			{

				//reconX[reconii]=zoom*(xx[fitnum*ii+jj]+left[ii]);
				//reconY[reconii]=zoom*(yy[fitnum*ii+jj]+top[ii]);

				//output section
				xbf[reconii]=xx[fitnum*ii+jj]+left[ii];
				topbf[reconii]=top[ii];
				leftbf[reconii]=left[ii];
				ybf[reconii]=yy[fitnum*ii+jj]+top[ii];
				bbf[reconii]=bbb[fitnum*ii+jj];
				tbf[reconii]=candz[ii];
				uncerxbf[reconii]=CRLBarray[(fitnum*2+1)*ii+jj];
				uncerybf[reconii]=CRLBarray[(fitnum*2+1)*ii+jj+1];
				nbf[reconii]=nn[fitnum*ii+jj];
				uncerbgbf[reconii]=CRLBarray[(fitnum*2+1)*ii+fittmp[ii]*2];
				covbf[reconii]=covariance[fitnum*ii+jj]/CRLBarray[(fitnum*2+1)*ii+jj]/CRLBarray[(fitnum*2+1)*ii+jj+1];
				NLLRbf[reconii]=div[ii];
				fitNbf[reconii]=(float) fittot[ii];


				//reconN[reconii]=nn[fitnum*ii+jj];
				//reconLAx[reconii]=zoom*CRLBarray[(fitnum*2+1)*ii+jj];
				////crlbx[reconii]=CRLBarray[(fitnum*2+1)*ii+jj];
				//reconLAy[reconii]=zoom*CRLBarray[(fitnum*2+1)*ii+jj+1];
				//reconcov[reconii]=covariance[fitnum*ii+jj]/CRLBarray[(fitnum*2+1)*ii+jj]/CRLBarray[(fitnum*2+1)*ii+jj+1];
				////cov[reconii]=covariance[fitnum*ii+jj]/CRLBarray[(fitnum*2+1)*ii+jj]/CRLBarray[(fitnum*2+1)*ii+jj+1];
				reconii++;
			}
		}
	}

	//repeat localization filtering.
	int tstart=0;
	int framereg=0;
	int modindex=0;
	float dis=-10;
	for(ii=0;ii<reconii;ii++)
	{
		framereg=(int) tbf[ii];
		if(xbf[ii]!=0)
		{
			for(jj=tstart;jj<reconii;jj++)
			{
				if(xbf[ii]==0) break;


				if(tbf[jj]!=framereg)
				{
					if(ii==jj)
						tstart=jj+1;
					break;
				}


				if((xbf[jj]!=0)&&(ii!=jj)
					&&(!((topbf[ii]==topbf[jj])&&(leftbf[ii]==leftbf[jj]))))
				{
					modindex=((uncerxbf[ii]+uncerybf[ii])<(uncerxbf[jj]+uncerybf[jj]))?jj:ii;
					dis=sqrt(pow(xbf[ii]-xbf[jj],2)+pow(ybf[ii]-ybf[jj],2));
					if (dis<=1*sqrt(pow(uncerxbf[modindex],2)+pow(uncerybf[modindex],2)))
					{
						xbf[modindex]=0;
						ybf[modindex]=0;
					}
				}
				

			}
		}
	}


	// counting procesure.
	int indextot=0;
	for(ii=0;ii<reconii;ii++){
		if(xbf[ii]!=0)
		{
			indextot=indextot+1;
			
		}
	}

	// allocation of reconstruction candidates


	float *reconX=(float*) malloc(indextot*sizeof(float));
	memset(reconX,0,indextot*sizeof(float));
	float *reconY=(float*) malloc(indextot*sizeof(float));
	memset(reconY,0,indextot*sizeof(float));
	float *reconN=(float*) malloc(indextot*sizeof(float));
	memset(reconN,0,indextot*sizeof(float));
	float *reconLAx=(float*) malloc(indextot*sizeof(float));
	memset(reconLAx,0,indextot*sizeof(float));
	float *reconLAy=(float*) malloc(indextot*sizeof(float));
	memset(reconLAy,0,indextot*sizeof(float));
	float *reconcov=(float*) malloc(indextot*sizeof(float));
	memset(reconcov,0,indextot*sizeof(float));

	// allocation of output candidates

	float * xout,* yout,* tout,* uncerxout,* unceryout;
	float *uncerbgout;
	float * covout,* NLLRout,*fitNout,*bout, *nout;
	
	outsize2[0]=indextot;
	outsize2[1]=1;
	plhs[1]= mxCreateNumericArray(2 ,outsize2, mxSINGLE_CLASS, mxREAL);
	xout=(float *)mxGetData(plhs[1]);

	plhs[2]= mxCreateNumericArray(2 ,outsize2, mxSINGLE_CLASS, mxREAL);
	yout=(float *)mxGetData(plhs[2]);

	plhs[3]= mxCreateNumericArray(2 ,outsize2, mxSINGLE_CLASS, mxREAL);
	tout=(float *)mxGetData(plhs[3]);

	plhs[4]= mxCreateNumericArray(2 ,outsize2, mxSINGLE_CLASS, mxREAL);
	uncerxout=(float *)mxGetData(plhs[4]);

	plhs[5]= mxCreateNumericArray(2 ,outsize2, mxSINGLE_CLASS, mxREAL);
	unceryout=(float *)mxGetData(plhs[5]);

	plhs[6]= mxCreateNumericArray(2 ,outsize2, mxSINGLE_CLASS, mxREAL);
	uncerbgout=(float *)mxGetData(plhs[6]);

	plhs[7]= mxCreateNumericArray(2 ,outsize2, mxSINGLE_CLASS, mxREAL);
	covout=(float *)mxGetData(plhs[7]);

	plhs[8]= mxCreateNumericArray(2 ,outsize2, mxSINGLE_CLASS, mxREAL);
	NLLRout=(float *)mxGetData(plhs[8]);

	plhs[9]= mxCreateNumericArray(2 ,outsize2, mxSINGLE_CLASS, mxREAL);
	fitNout=(float*)mxGetData(plhs[9]);

	plhs[10]= mxCreateNumericArray(2 ,outsize2, mxSINGLE_CLASS, mxREAL);
	bout=(float *)mxGetData(plhs[10]);

	plhs[11]=mxCreateNumericArray(2 ,outsize2, mxSINGLE_CLASS, mxREAL);
	nout=(float *)mxGetData(plhs[11]);


	//record coordinates for output and reconstruction
	int recdx=0;
	for(ii=0;ii<reconii;ii++)
	{
		if(xbf[ii]!=0)
		{ //record
			//reconstruction candidates
			reconX[recdx]=xbf[ii]*zoom;
			reconY[recdx]=ybf[ii]*zoom;
			reconN[recdx]=nbf[ii];
			reconLAx[recdx]=uncerxbf[ii]*zoom;
			reconLAy[recdx]=uncerybf[ii]*zoom;
			reconcov[recdx]=covbf[ii];

			//out put candidates
			xout[recdx]=xbf[ii];
			yout[recdx]=ybf[ii];
			tout[recdx]=tbf[ii];
			uncerxout[recdx]=uncerxbf[ii];
			unceryout[recdx]=uncerybf[ii];
			uncerbgout[recdx]=uncerbgbf[ii];
			covout[recdx]=covbf[ii];
			NLLRout[recdx]=NLLRbf[ii];
			fitNout[recdx]=fitNbf[ii];
			bout[recdx]=bbf[ii];
			nout[recdx]=nbf[ii];
			recdx++;
		}
	}



	printTiming("t3",t3);
	QueryPerformanceCounter(&t4);

	//reconstruction! 

	outsize[0]=(int) floor(data_xsz*zoom);
	outsize[1]=(int) floor(data_ysz*zoom);
	plhs[0]= mxCreateNumericArray(2, outsize, mxSINGLE_CLASS, mxREAL);
	SuperResim=(float *)mxGetData(plhs[0]);
	int xsz=(int) outsize[0];
	int ysz=(int) outsize[1];

	/*outsize2[0]=20;
	outsize2[1]=20;
	outsize2[2]=reconii;
	plhs[1]= mxCreateNumericArray(Adims ,outsize2, mxSINGLE_CLASS, mxREAL);
	out=(float *)mxGetData(plhs[1]);*/ 


	GPUgenerateblobs(indextot, xsz, ysz,reconX,reconY, reconN, reconLAx,reconLAy,reconcov,SuperResim);
	printTiming("t4", t4);

	free(CRLBarray);
	free(covariance);
	free(reconX);
	free(reconY);
	free(reconN);
	free(reconLAx);
	free(reconLAy);
	free(reconcov);
	free(left);
	free(top);
	free(in_sub);
	free(bbb);
	free(div);
	free(xx);
	free(yy);
	free(nn);
	free(fittot);
	free(fittmp);

	free(xbf);
	free(ybf);
	free(tbf);
	free(nbf);
	free(bbf);
	free(uncerxbf);
	free(uncerybf);
	free(uncerbgbf);
	free(covbf);
	free(NLLRbf);
	free(fitNbf);

}