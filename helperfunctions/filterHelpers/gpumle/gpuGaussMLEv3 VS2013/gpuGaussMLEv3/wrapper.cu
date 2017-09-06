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
extern "C" void kernel_MLERatio_wrapper(dim3 dimGrid, dim3 dimBlock, cudaStream_t stream, float *d_data, float *PSFSigma, int szXY, int szZ, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits){

	if (szZ == 1)
		kernel_MLERatio << <dimGrid, dimBlock, 0, stream >> >(d_data, PSFSigma[0], szXY, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
	/*else if (szZ > 1)
	kernel_MLEFit3D << <dimGrid, dimBlock, 0, stream >> >(d_data, PSFSigma[0], PSFSigma[1], szXY, szZ, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);*/

}

//*******************************************************************************************
extern "C" void kernel_MLEFit_wrapper(dim3 dimGrid, dim3 dimBlock, cudaStream_t stream, float *d_data, float *PSFSigma, int szXY, int szZ, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits){
	/*!
	*  \brief Basic maximum likelihood estimator fit based kernel
	*  \param dimGrid number of blocks
	*  \param dimBlock number of threads per block
	*  \param d_data an array of subregions to be processed copied into video memory
	*  \param PSFSigma the sigma value to use for the point spread function
	*  \param sz nxn size of the subregion
	*  \param iterations maximum allowed iterations before aborting fitting
	*  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
	*  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates
	*  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
	*  \param Nfits number of subregions to fit
	*/
	
	if (szZ == 1)
		kernel_MLEFit << <dimGrid, dimBlock, 0, stream >> >(d_data, PSFSigma[0], szXY, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
	else if(szZ > 1)
		kernel_MLEFit3D << <dimGrid, dimBlock, 0, stream >> >(d_data, PSFSigma[0], PSFSigma[1], szXY, szZ, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
		
}


extern "C" void kernel_MLEFitSigma_wrapper(dim3 dimGrid, dim3 dimBlock, cudaStream_t stream, float *d_data, float *PSFSigma, int szXY, int szZ, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits){

	if (szZ == 1)
		kernel_MLEFitSigma << <dimGrid, dimBlock, 0, stream >> >(d_data, PSFSigma[0], szXY, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
	else if (szZ > 1)
		kernel_MLEFit3DSigma << <dimGrid, dimBlock, 0, stream >> >(d_data, PSFSigma[0], PSFSigma[1], szXY, szZ, iterations,
			d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);

}

extern "C" void kernel_MLEFitSigmaXY_wrapper(dim3 dimGrid, dim3 dimBlock, cudaStream_t stream, float *d_data, float *PSFSigma, int szXY, int szZ, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits){

	//if (szZ == 1)
	kernel_MLEFitSigmaXX << <dimGrid, dimBlock, 0, stream >> >(d_data, PSFSigma[0], szXY, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
}