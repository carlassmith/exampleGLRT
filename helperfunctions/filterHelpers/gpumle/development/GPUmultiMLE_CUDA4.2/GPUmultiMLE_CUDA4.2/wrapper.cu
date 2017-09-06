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
extern "C" void kernel_example_wrapper(dim3 dimGrid, dim3 dimBlock, int *c, const int *a, const int *b)
{
	/*! 
	 * \brief Example kernel to test build settings.
	 */
	kernel_example<<<dimGrid,dimBlock>>>(c,a,b);
}