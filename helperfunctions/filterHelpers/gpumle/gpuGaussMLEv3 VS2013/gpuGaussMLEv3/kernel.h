/*!
 * \file kernel.h
 * \author <your name>
 * \date April 20, 2012
 * \brief  Put prototypes of all Cuda Kernels here.
 */

#ifndef KERNEL_H
#define KERNEL_H
#include <cuda_runtime.h>

__global__ void kernel_MLEFit(const float *d_data, float PSFSigma, int sz, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits);

__global__ void kernel_MLEFit3D(const float *d_data, float PSFSigmax, float PSFSigmaz , int szXY, int szZ, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits);

__global__ void kernel_MLEFitSigma(float *d_data, float PSFSigma, int sz, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits);

__global__ void kernel_MLEFit3DSigma(const float *d_data, float PSFSigmax, float PSFSigmaz, int szXY, int szZ, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits);

__global__ void kernel_MLERatio(const float *d_data, float PSFSigma, int sz, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits);

__global__ void kernel_MLEFitSigmaXX(const float *d_data, float PSFSigma, int sz, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits);

#endif