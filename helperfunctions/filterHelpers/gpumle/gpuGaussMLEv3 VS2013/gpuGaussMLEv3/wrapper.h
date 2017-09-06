

extern "C" void kernel_MLEFit_wrapper(dim3 dimGrid, dim3 dimBlock, cudaStream_t stream, float *d_data, float *PSFSigma, int szXY, int szZ, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits);


extern "C" void kernel_MLEFitSigma_wrapper(dim3 dimGrid, dim3 dimBlock, cudaStream_t stream, float *d_data, float *PSFSigma, int szXY, int szZ, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits);

extern "C" void kernel_MLEFitSigmaXY_wrapper(dim3 dimGrid, dim3 dimBlock, cudaStream_t stream, float *d_data, float *PSFSigma, int szXY, int szZ, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits);

extern "C" void kernel_MLERatio_wrapper(dim3 dimGrid, dim3 dimBlock, cudaStream_t stream, float *d_data, float *PSFSigma, int szXY, int szZ, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits);