#include <windows.h>
#pragma comment(lib, "kernel32.lib")

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
#include <cuda_runtime.h>
#include "definitions.h"

#ifndef max
//! not defined in the C standard used by visual studio
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
//! not defined in the C standard used by visual studio
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

//******************************************************
// extern C declarations of the kernels from wrapper.cu.  Yes you can put this in a
// header if you want.
extern "C" void kernel_example_wrapper(dim3 dimGrid, dim3 dimBlock, int* c, const int* a, const int*b);

//*******************************************************************************************
void cudaerror() {
/*!
 *  \brief A simple function to dump Cuda errors and abort the run.
 */
    cudaError_t errornum;
    const char *str=0;
    if (errornum=cudaGetLastError()) {
        str=cudaGetErrorString(errornum);
		mexErrMsgIdAndTxt("CudaTemplate:CUDA","%s\nYou should clear this function in MATLAB for proper operation.\n", str);
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
	//declare all vars
	int *array_a=0, *array_b=0; 
	int *d_a=0, *d_b=0, *d_c=0;
	const mwSize *array_a_size=0,*array_b_size=0;
	size_t Ndim_a=0, Ndim_b=0;

	//check for required inputs, correct types, and dimensions
	if (nrhs != 2)
		mexErrMsgIdAndTxt("CUDAtemplate:WrongNumberInputs","Input must include vectors a,b");

	if (mxGetClassID(prhs[0])!=mxINT32_CLASS)
		mexErrMsgIdAndTxt("CUDAtemplate:WrongType","Inputs must be 32 bit ints");
	
	Ndim_a = mxGetNumberOfDimensions(prhs[0]);
	Ndim_b = mxGetNumberOfDimensions(prhs[1]);
	array_a_size = mxGetDimensions(prhs[0]);
	array_b_size = mxGetDimensions(prhs[1]);

	//1D vectors still return 2D
	if (Ndim_a > 2 || Ndim_b > 2 || array_a_size[0] > 1 || array_b_size[0] > 1)
		mexErrMsgIdAndTxt("CUDAtemplate:WrongDimensions","a and b should be 1D vectors");

	if (array_a_size[1] != array_b_size[1])
		mexErrMsgIdAndTxt("CUDAtemplate:VectorMismatch","size of a != b");
	//retrieve all inputs
	array_a = (int*) mxGetData(prhs[0]);
	array_b = (int*) mxGetData(prhs[1]);

	//validate input values(this section better not be blank!)

	//allocate memory and shove data on the GPU
	dim3 dimBlock((unsigned int) array_a_size[1]);
    dim3 dimGrid(1);
	const size_t bytes = array_a_size[1]*sizeof(float);
	cudaMalloc((void**)&d_a, bytes);
	cudaMalloc((void**)&d_b, bytes);
	cudaMalloc((void**)&d_c, bytes);
	cudaMemcpy(d_a, array_a, bytes, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, array_b, bytes, cudaMemcpyHostToDevice);

	//kernel calls
	kernel_example_wrapper(dimGrid,dimBlock,d_c,d_a,d_b);
	cudaerror();

	//allocate results memory and copy results back to matlab
	plhs[0] = mxCreateNumericMatrix(array_a_size[0],array_a_size[1],mxINT32_CLASS, mxREAL);
	cudaMemcpy((int*)mxGetData(plhs[0]), d_c, bytes, cudaMemcpyDeviceToHost);

	//clean up!
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);
	cudaThreadExit(); //release context so future cudaSetDevice calls work

	return;
 }