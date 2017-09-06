/* sample usage
  a = int32(randi(1000,1000,1))';
  b = int32(randi(1000,1000,1))';
  c = mexfilename(a,b);
*/

//#include <windows.h>
//#pragma comment(lib, "kernel32.lib")

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
#include <timer.h>
#include "kernel.h"

// include CUDA
#include <cuda_runtime.h>

// includes, project
#include <helper_cuda.h>
#include <helper_functions.h>  // helper for shared that are common to CUDA Samples

#include "definitions.h"
#include "wrapper.h"

#ifndef max
//! not defined in the C standard used by visual studio
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
//! not defined in the C standard used by visual studio
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

//*******************************************************************************************
void CUDAERRROR(const char *instr) {
	/*!
	*  \brief A simple function to dump Cuda errors and abort the run.
	*  \param instr unused string pointer
	*/
	cudaError_t errornum;
	const char *str;
	if ((errornum = cudaGetLastError())) {
		str = cudaGetErrorString(errornum);
		mexPrintf("%s:%s\n", instr, str);
		mexPrintf("You should clear this function in MATLAB for proper operation.\n", str);
	}
}



//struct numberOfVarFitType lookupNumberOfVarFitType[] = { { , , }, // (bg),(bg,I)
//{  }, // (bg),(x,y,bg,I)
//{ NV_PSH0, NV_PSH1, Q_PS }, // (bg),(x,y,bg,I,Sigma)
//{ NV_PH03, NV_PH13, Q_P3 }, // (bg), (x,y,bg,I,z)
//{ 1, 6, 6 }, // (bg),(x,y,bg,I,Sx,Sy)	
//{ NV_PH05, NV_PH15, Q_P5 }//(bg)(x,y,z,bg,I,Sx,Sz)
//};

struct numberOfVarFitType{
	size_t NParametersToFitH0 = 0;
	size_t NParametersToFitH1 = 0;
	size_t NQParametersOfFit = 0;
};

//const numberOfVarFitType lookupNumberOfVarFitType[] = { { }, // (bg),(bg,I)
//{ NV_PH0, NV_PH1, Q_P }, // (bg),(x,y,bg,I)
//{ NV_PSH0, NV_PSH1, Q_PS }, // (bg),(x,y,bg,I,Sigma)
//{ NV_PH03, NV_PH13, Q_P3 }, // (bg), (x,y,bg,I,z)
//{ 1, 6, 6 }, // (bg),(x,y,bg,I,Sx,Sy)
//{  }//(bg)(x,y,z,bg,I,Sx,Sz)
//};


numberOfVarFitType getLookupNumberOfVarFitTypeStruct(int fitType,size_t ZetSlices){

	
	numberOfVarFitType newStruct;

	switch (fitType){
	case 0:
		newStruct.NParametersToFitH0 = NV_RH0;
		newStruct.NParametersToFitH1 = NV_RH1;
		newStruct.NQParametersOfFit = Q_R;
		break;
	case 1:
		newStruct.NParametersToFitH0 = NV_PH0;
		newStruct.NParametersToFitH1 = NV_PH1;
		newStruct.NQParametersOfFit = Q_P;
		break;
	case 2:
		newStruct.NParametersToFitH0 = NV_PSH0;
		newStruct.NParametersToFitH1 = NV_PSH1;
		newStruct.NQParametersOfFit = Q_PS;
		break;
	case 3:
		newStruct.NParametersToFitH0 = (NV_PH03 +  ZetSlices - 1);
		newStruct.NParametersToFitH1 = (NV_PH13 + ZetSlices - 1);
		newStruct.NQParametersOfFit = Q_P3;
		break;
	case 4:
		newStruct.NParametersToFitH0 = 1;
		newStruct.NParametersToFitH1 = 6;
		newStruct.NQParametersOfFit = 6;
		break;
	case 5:
		newStruct.NParametersToFitH0 = NV_PH05;
		newStruct.NParametersToFitH1 = NV_PH15;
		newStruct.NQParametersOfFit = Q_P5;
		break;
	default:
		newStruct.NParametersToFitH0 = 0;
		newStruct.NParametersToFitH1 = 0;
		newStruct.NQParametersOfFit = 0;
	break;
	}
	return newStruct;
}




void getDeviceFitProberties(TGPUplan* plan, int devCount, double *gpuDist, size_t szXY, size_t szZ, bool showProp){

	bool limitFits = ((gpuDist != NULL) && gpuDist[0] > 0);
	int summedFits = 0;
	for (int i = 0; i < devCount; ++i){
		cudaDeviceProp devProp;
		checkCudaErrors(cudaGetDeviceProperties(&devProp, i));

		plan[i].MAXMEMS = (devProp.sharedMemPerBlock) / sizeof(float);
		int BlockSize = (int)floor((float)plan[i].MAXMEMS / (float)szXY / (float)szXY / (float)szZ / (float)devProp.warpSize) * devProp.warpSize;
		plan[i].BlockSize = min(devProp.maxThreadsPerBlock, max(4, BlockSize));

		plan[i].IMBIGS = (int)floor(((float) sizeof(float)* plan[i].MAXMEMS / devProp.warpSize));

		plan[i].maxGrid = (int)min(floor((float)devProp.totalGlobalMem / ((float)devProp.sharedMemPerBlock) / (2.0f)), devProp.maxGridSize[0]);
		if (limitFits)
			plan[i].maxNumberOfFits = (int)min((double)plan[i].maxGrid *plan[i].BlockSize, gpuDist[i]);
		else
			plan[i].maxNumberOfFits = plan[i].maxGrid *plan[i].BlockSize;
		summedFits += (int)plan[i].maxNumberOfFits;

		plan[i].MEMS = BlockSize*szXY*szXY*szZ;

		//this variable is double
		plan[i].dimBlock = plan[i].BlockSize;
		plan[i].major = devProp.major;
		strcpy(plan[i].name, devProp.name);

		if (showProp){
			mexPrintf("currentImSize: %dx%dx%d\n", (int)szXY, (int)szXY, (int)szZ);
			mexPrintf("plan[%d].name %s\n", i, plan[i].name);
			mexPrintf("plan[%d].numberOfThreadsPerBlock %d\n\n", i, plan[i].BlockSize);
			mexPrintf("plan[%d].maxNumberOfFitsInGlobalMem %d\n", i, plan[i].maxNumberOfFits);
			//mexPrintf("plan[%d].maxAvailableMemPerBlock %d\n", i, plan[i].MAXMEMS);
			mexPrintf("plan[%d].maxPixPerThread %d\n", i, plan[i].IMBIGS);
			mexPrintf("plan[%d].MaxGridSize %d\n", i, plan[i].maxGrid);
			mexPrintf("plan[%d].MaxUsedMemPerBlock %d\n", i, plan[i].MEMS);
			mexPrintf("plan[%d].computeCapabilityMajor %d\n", i, plan[i].major);
		}
	}
	if (summedFits == 0) mexErrMsgIdAndTxt("gpuDist:Invalid parameter", "Fit limit is equal 0. Check gpuDist");

}

void createAllDeviceStreams(TGPUplan* plan, int devCount,bool show){
	//Create Streams to devices
	for (int i = 0; i < devCount; i++)
	{
			checkCudaErrors(cudaSetDevice(i));
			if (show) mexPrintf("Create device stream (%d)\n", i);
			checkCudaErrors(cudaStreamCreate(&plan[i].stream));
	
	}
}

void destroyAllDeviceStreams(TGPUplan* plan, int devCount,bool show){

	//Shut down this GPU
	for (int i = 0; i < devCount; i++)
	{
			//Shut down this GPU
			checkCudaErrors(cudaSetDevice(i));
			if (show) mexPrintf("Destroy device stream (%d)\n", i);
			checkCudaErrors(cudaStreamDestroy(plan[i].stream));
			cudaDeviceReset();
	}
}
//*******************************************************************************************
void mexFunction(int nlhs, mxArray *plhs[], int	nrhs, const	mxArray	*prhs[]) {
	/*!
		*  \brief Entry point in the code for Matlab.  Equivalent to main().
		*  \param nlhs number of left hand mxArrays to return
		*  \param plhs array of pointers to the output mxArrays
		*  \param nrhs number of input mxArrays
		*  \param prhs array of pointers to the input mxArrays.
		*/

	size_t  Nfitraw, szXY, szZ;
	bool show;

	// Get dimentions from Matlab
	const mwSize *datasize = mxGetDimensions(prhs[0]);
	size_t Ndim = mxGetNumberOfDimensions(prhs[0]);

	const mwSize * sigmaSize = mxGetDimensions(prhs[1]);
	size_t NdimSigma = mxGetNumberOfDimensions(prhs[1]);
	
	//Get data
	float *data = (float *)mxGetData(prhs[0]);
	float *PSFSigma = (float*)mxGetData(prhs[1]);
	int iterations = (int)mxGetScalar(prhs[2]);
	size_t fittype = (int)mxGetScalar(prhs[3]);

	// Get Device count
	TGPUplan* plan;
	int devCount;
	cudaGetDeviceCount(&devCount);
	if (devCount == 0)
		mexPrintf("No CUDA enabled devices detected.");

	//Get data size based on fittype
	switch (Ndim){
	case 2:
		Nfitraw = 1;
		szXY = datasize[0];
		szZ = 1;
		break;
	case 3:
		szXY = datasize[0];
		switch (fittype)
		{
		case 3: case 5:
			Nfitraw = 1;
			szZ = datasize[2];
			break;
		case 0: case 1: case 2: case 4:
			Nfitraw = datasize[2];
			szZ = 1;
			break;
		default:
			mexErrMsgIdAndTxt("gpuGaussMLE:Invalid parameter", "This fittype not implemented");
		}
		break;
	case 4:
		szXY = datasize[0];
		szZ = datasize[2];
		Nfitraw = datasize[3];
		break;
	default:
		mexErrMsgIdAndTxt("gpuGaussMLE:WrongDimensions","Data dimentions not supported");
	}

	numberOfVarFitType lookupNumberOfVarFitType = getLookupNumberOfVarFitTypeStruct(fittype, szZ);

	double *gpuDist = NULL;
	size_t gpuDistNdim;
	const mwSize *gpuDistSize;

	//Check input
	if (!mxIsSingle(prhs[0]) || !mxIsSingle(prhs[1]))
		mexErrMsgTxt("Data and PSFSigma must be single arrays.");

	if (Ndim > 2 && (NdimSigma < 1 || NdimSigma > 2 || sigmaSize[0] != 1)) //min 2 dim
		mexErrMsgIdAndTxt("gpuGaussMLE:WrongDimensions", "Sigma should be a 1xn dimensional array.\n");

	if (szZ > 1 && sigmaSize[1] != 2) //min 2 dim
		mexErrMsgIdAndTxt("gpuGaussMLE:WrongDimensions", "Sigma should be a 1x2 vector for this fittype.\n");
	if (szZ == 1 && sigmaSize[1] != 1) //min 2 dim
		mexErrMsgIdAndTxt("gpuGaussMLE:WrongDimensions", "Sigma should be a scalar for this fittype.\n");



	if (nrhs >= 5)
		show = ((double) mxGetScalar(prhs[4]) > 0);
	else 
		show = false;

	switch (nrhs)
	{
		case 4: case 5:
			if(show) mexPrintf("Equal distribution over CUDA enabled GPUs.\n");  //Run function with gpuDist = -1 to list all CUDA enabled GPUs and there properties
			break;
		case 6:
			gpuDistSize = mxGetDimensions(prhs[5]);
			gpuDistNdim = mxGetNumberOfDimensions(prhs[5]);
			gpuDist = (double*)mxGetData(prhs[5]);
			if (gpuDist[0] >= 0 && ((gpuDistNdim > 1 && (gpuDistNdim < 1 || gpuDistNdim > 2 || gpuDistSize[0] != 1)) || gpuDistSize[1] != devCount)){//min 2 dim
				if (gpuDistSize[1] != devCount)
					mexPrintf("\n\nThere are %d CUDA enabled GPUs detected.To list all CUDA enabled GPUs try gpuDist = -1. \n", devCount);
			
				mexErrMsgIdAndTxt("gpuDist:WrongDimensions", "gpuDist should be a 1xn dimensional array with maximum number of fits per GPU, where n is the number of CUDA enabled devices.\n");
			}
			
			break;
		default:
			mexErrMsgIdAndTxt("gpuGaussMLE:WrongDimensions", "Wrong number input parameters\n");
			break;
	}

	// Setup output Matlab variables
	plhs[0] = mxCreateNumericMatrix((lookupNumberOfVarFitType.NParametersToFitH0 + lookupNumberOfVarFitType.NParametersToFitH1), Nfitraw, mxSINGLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(lookupNumberOfVarFitType.NParametersToFitH1, Nfitraw, mxSINGLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericMatrix(lookupNumberOfVarFitType.NQParametersOfFit, Nfitraw, mxSINGLE_CLASS, mxREAL);

	//create device variable for data and copy to device
	if (show || ((gpuDist != NULL) && gpuDist[0] < 0)) printf("CUDA-capable device count: %i\n", devCount);

	plan = (TGPUplan *)malloc(sizeof(TGPUplan)*devCount);

	getDeviceFitProberties(plan, devCount, gpuDist, szXY, szZ, (show || ((gpuDist != NULL) && gpuDist[0] < 0)));
	for (int i = 0; i < devCount; i++)
	{
		if (plan[i].major < 2){ //TODO && plan[i].maxNumberOfFits > 0
			char err[400];
			sprintf(err, "GPU device: %s is too old.\nUpdate your GPU device (compute capability > 2.0).", plan[i].name);
			mexErrMsgTxt(err);
			return;
		}
	}
	//mexPrintf("2.2f\n", gpuDist[0]);
	// Check if we stop here
	if ((gpuDist != NULL) && gpuDist[0] < 0){
		mexPrintf("Listed all CUDA enabled GPUs. Done.\n");
		return;
	}

	//create GPU streams to CUDA enabled devices
	createAllDeviceStreams(plan, devCount,show);

	size_t placedFits = 0;
	mexEvalString("drawnow");

	placedFits = 0;
	int currentFit = 0;
	int copiedFits = 0;
	int devicesUsed = devCount;
	int circleFits = 0;
	while (placedFits < Nfitraw){
		circleFits = 0;
		//Create streams for issuing GPU command asynchronously and allocate memory 
		mexPrintf("Allocating host and device memory...\n");
		mexEvalString("drawnow");
		for (int i = 0; i < devCount; i++)
		{
			plan[i].Nfits = (int)min(Nfitraw - placedFits, plan[i].maxNumberOfFits);
			plan[i].dimGrid = (int)min(ceil((float)plan[i].Nfits / (float)plan[i].BlockSize), plan[i].maxGrid);
			plan[i].Nfits = min(plan[i].Nfits, plan[i].BlockSize*plan[i].dimGrid);

			placedFits += plan[i].Nfits;
			circleFits += plan[i].Nfits;
			mexPrintf("plan[%d].Nfits %d (blockDim %d, gridDim %d) on %s\n", i, plan[i].Nfits, plan[i].dimBlock, plan[i].dimGrid,plan[i].name);
			
			//Allocate memory
			checkCudaErrors(cudaSetDevice(i));
			checkCudaErrors(cudaMalloc((void**)&plan[i].d_Data, szXY*szXY*szZ*plan[i].Nfits*sizeof(float)));

			checkCudaErrors(cudaMalloc((void**)&plan[i].d_Parameters,
				(lookupNumberOfVarFitType.NParametersToFitH0 + lookupNumberOfVarFitType.NParametersToFitH1)*plan[i].Nfits*sizeof(float)));
			checkCudaErrors(cudaMalloc((void**)&plan[i].d_CRLBs, lookupNumberOfVarFitType.NParametersToFitH1*plan[i].Nfits*sizeof(float)));
			checkCudaErrors(cudaMalloc((void**)&plan[i].d_LogLikelihood, lookupNumberOfVarFitType.NQParametersOfFit*plan[i].Nfits*sizeof(float)));

			checkCudaErrors(cudaMallocHost((void **)&plan[i].h_Parameters,
				(lookupNumberOfVarFitType.NParametersToFitH0 + lookupNumberOfVarFitType.NParametersToFitH1)*plan[i].Nfits*sizeof(float)));
			checkCudaErrors(cudaMallocHost((void**)&plan[i].h_CRLBs, lookupNumberOfVarFitType.NParametersToFitH1*plan[i].Nfits*sizeof(float)));
			checkCudaErrors(cudaMallocHost((void**)&plan[i].h_LogLikelihood, lookupNumberOfVarFitType.NQParametersOfFit*plan[i].Nfits*sizeof(float)));
			
			//cudaMemsetAsync
			//checkCudaErrors(cudaMallocHost((void**)&plan[i].h_Data, szXY*szXY*szZ*plan[i].Nfits*sizeof(float)));
			if ((Nfitraw - placedFits) == 0){
				devicesUsed = i + 1;
				break;
			}
		}

		//Start timing and compute on GPU(s)
		if (placedFits < Nfitraw)
			mexPrintf("Placed %d of %d fits on the GPU %d devices\n", placedFits, Nfitraw, devicesUsed);
		else
			mexPrintf("All fits are distributed over the GPU %d devices\n", devicesUsed);
		mexEvalString("drawnow");
		StartTimer();
		
		//Copy data to GPU, launch the kernel and copy data back. All asynchronously
		for (int i = 0; i < devicesUsed; i++)
		{
			//Set device
			checkCudaErrors(cudaSetDevice(i));
			// saves one ms not to first put it in the "optimized" cudaMallocHost.
			//checkCudaErrors(cudaMemcpyAsync(plan[i].h_Data, data + szXY*szXY*szZ*currentFit, szXY*szXY*szZ*plan[i].Nfits*sizeof(float), cudaMemcpyHostToHost, plan[i].stream));
			//checkCudaErrors(cudaMemcpyAsync(plan[i].d_Data, plan[i].h_Data, szXY*szXY*szZ*plan[i].Nfits * sizeof(float), cudaMemcpyHostToDevice, plan[i].stream));
			
			checkCudaErrors(cudaMemcpyAsync(plan[i].d_Data, data + szXY*szXY*szZ*currentFit, szXY*szXY*szZ*plan[i].Nfits * sizeof(float), cudaMemcpyHostToDevice, plan[i].stream));
			currentFit += (int)plan[i].Nfits;
			
			dim3 dimB(plan[i].dimBlock);
			dim3 dimG(plan[i].dimGrid);

			//Perform GPU computations
			switch (fittype)
			{
				case 0:
					kernel_MLERatio_wrapper(dimG, dimB, plan[i].stream, plan[i].d_Data, PSFSigma, (int)szXY, (int)szZ, iterations,
						plan[i].d_Parameters, plan[i].d_CRLBs, plan[i].d_LogLikelihood, (int)plan[i].Nfits);
					break;
				case 1: case 3:
					kernel_MLEFit_wrapper(dimG, dimB, plan[i].stream, plan[i].d_Data, PSFSigma, (int)szXY, (int)szZ, iterations,
					plan[i].d_Parameters, plan[i].d_CRLBs, plan[i].d_LogLikelihood, (int)plan[i].Nfits);
					break;
				case 2: case 5:
					kernel_MLEFitSigma_wrapper(dimG, dimB, plan[i].stream, plan[i].d_Data, PSFSigma, (int) szXY, (int) szZ, iterations,
						plan[i].d_Parameters, plan[i].d_CRLBs, plan[i].d_LogLikelihood, (int)plan[i].Nfits);
					break;
				case 4:
						kernel_MLEFitSigmaXY_wrapper(dimG, dimB, plan[i].stream, plan[i].d_Data, PSFSigma, (int)szXY, (int)szZ, iterations,
						plan[i].d_Parameters, plan[i].d_CRLBs, plan[i].d_LogLikelihood, (int)plan[i].Nfits);
					break;

			default:
				break;
			}
			/*mexPrintf("%d", (int)szZ);
			mexPrintf("%d", (int)szXY);
			mexPrintf("%d", (int)Ndim);*/
			CUDAERRROR("reduceKernel() execution failed.\n");
			
			//Read back GPU results
			checkCudaErrors(cudaMemcpyAsync(plan[i].h_Parameters, plan[i].d_Parameters, 
				(lookupNumberOfVarFitType.NParametersToFitH0 + lookupNumberOfVarFitType.NParametersToFitH1)*plan[i].Nfits*sizeof(float), cudaMemcpyDeviceToHost, plan[i].stream));
			checkCudaErrors(cudaMemcpyAsync(plan[i].h_CRLBs, plan[i].d_CRLBs, 
				lookupNumberOfVarFitType.NParametersToFitH1*plan[i].Nfits*sizeof(float), cudaMemcpyDeviceToHost, plan[i].stream));
			checkCudaErrors(cudaMemcpyAsync(plan[i].h_LogLikelihood, plan[i].d_LogLikelihood, lookupNumberOfVarFitType.NQParametersOfFit*plan[i].Nfits*sizeof(float), cudaMemcpyDeviceToHost, plan[i].stream));
			
		}
		mexPrintf("Cleanup...\n\n");
		mexEvalString("drawnow");

		//Process GPU results
		for (int i = 0; i < devicesUsed; i++)
		{
			//Set device
			checkCudaErrors(cudaSetDevice(i));

			//Wait for all operations to finish
			cudaStreamSynchronize(plan[i].stream);

			//Merge Host arrays		
			memcpy((float *)mxGetData(plhs[0]) + copiedFits* (lookupNumberOfVarFitType.NParametersToFitH1 + lookupNumberOfVarFitType.NParametersToFitH0),
				plan[i].h_Parameters, plan[i].Nfits* (lookupNumberOfVarFitType.NParametersToFitH1 + lookupNumberOfVarFitType.NParametersToFitH0)*sizeof(float));
			memcpy((float *)mxGetData(plhs[1]) + copiedFits*lookupNumberOfVarFitType.NParametersToFitH1, plan[i].h_CRLBs, plan[i].Nfits*lookupNumberOfVarFitType.NParametersToFitH1*sizeof(float));
			memcpy((float *)mxGetData(plhs[2]) + copiedFits*lookupNumberOfVarFitType.NQParametersOfFit, plan[i].h_LogLikelihood, plan[i].Nfits*lookupNumberOfVarFitType.NQParametersOfFit*sizeof(float));
			copiedFits += plan[i].Nfits;
			

			// Cleanup
			checkCudaErrors(cudaFree(plan[i].d_Data));
			checkCudaErrors(cudaFree(plan[i].d_CRLBs));
			checkCudaErrors(cudaFree(plan[i].d_LogLikelihood));
			checkCudaErrors(cudaFree(plan[i].d_Parameters));

			//checkCudaErrors(cudaFreeHost(plan[i].h_Data));
			checkCudaErrors(cudaFreeHost(plan[i].h_CRLBs));
			checkCudaErrors(cudaFreeHost(plan[i].h_LogLikelihood));
			checkCudaErrors(cudaFreeHost(plan[i].h_Parameters));	
			
		}
		mexPrintf("  GPU Processing time %f (ms) for %d fits\n\n", GetTimer(), circleFits);
		mexEvalString("drawnow");
	}

	destroyAllDeviceStreams(plan, devCount, show);
	
	return;
}
