/*!
 * \file GPUgaussLib.h
 * \author Keith Lidke
 * \date January 10, 2010
 * \brief Prototypes for all the Cuda helper functions.
 */

// This code provides a set of functions that can be called from inside
// NVIDIA CUDA Kernels.
__device__ void avgDerivative(int, float *, float * aX, float * aY);

__device__ float kernel_IntGauss1D(const int ii, const float x, const float sigma);

__device__ float kernel_alpha(const float z, const float Ax, const float Bx, const float d);

__device__ float kernel_dalphadz(const float z, const float Ax, const float Bx, const float d);

__device__ float kernel_d2alphadz2(const float z, const float Ax, const float Bx, const float d);


__device__ void kernel_CenterofMass2D(const int sz, const float *data, float *x, float *y, float aX, float aY);

__device__ void kernel_GaussFMaxMin2D(const int sz, const float sigma, const float * data, float *MaxN, float *MinBG);