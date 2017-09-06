/*!
 * \file wrapper.cu
 * \author 
 * \date April 20, 2012
 * \brief Wrap the Cuda kernel calls as standard external C functions.  This allows the kernels to be
 * called without doing anything special in the C code and simplifies building the code.
 */

#include "definitions.h"
#include "kernel.h"


//*******************************************************************************************
extern void kernel_guassiansampleblobs_wrapper(dim3 dimGrid, dim3 dimBlock, int iiK,int BlockSize, int sz, float *d_xarray,float *d_yarray,float *d_Narray, float *d_xsigma,float *d_ysigma,float *d_covariance,float *d_im,float *d_xl,float *d_yt  ) 
{
	/*! 
	 * \brief Example kernel to test build settings.
	 */
	kernel_guassiansampleblobs<<<dimGrid, dimBlock>>>(iiK,BlockSize,sz, d_xarray,d_yarray,d_Narray, d_xsigma,d_ysigma,d_covariance,d_im,d_xl,d_yt);
}


extern void kernel_guassianintegrateblobs_wrapper(dim3 dimGrid, dim3 dimBlock, int iiK,int BlockSize, int sz, float *d_xarray,float *d_yarray,float *d_Narray, float *d_xsigma,float *d_ysigma,float *d_covariance,float *d_im,float *d_xl,float *d_yt  ) 
{
	/*! 
	 * \brief Example kernel to test build settings.
	 */
	kernel_guassianintegrateblobs<<<dimGrid, dimBlock>>>(iiK,BlockSize,sz, d_xarray,d_yarray,d_Narray, d_xsigma,d_ysigma,d_covariance,d_im,d_xl,d_yt);
}