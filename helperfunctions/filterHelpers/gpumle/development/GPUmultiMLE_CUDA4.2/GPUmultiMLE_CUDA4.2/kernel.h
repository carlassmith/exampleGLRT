/*!
 * \file kernel.h
 * \author <your name>
 * \date April 20, 2012
 * \brief  Put prototypes of all Cuda Kernels here.
 */

#ifndef KERNEL_H
#define KERNEL_H

__global__ void kernel_example(int *c, const int *a, const int *b);

#endif