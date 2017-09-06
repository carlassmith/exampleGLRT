/*!
 * \file kernel.h
 * \author <your name>
 * \date April 20, 2012
 * \brief  Put prototypes of all Cuda Kernels here.
 */

#ifndef KERNEL_H
#define KERNEL_H
#define NK 128 //number of blocks to run in each kernel
// Thread block size
#define BSZ 64
#define MEM 70
#define IMSZ 11
#define IMSZBIG 21
#define imMEM 4000


__global__ void kernel_guassiansampleblobs(int,int,int, float*,float*,float*, float*,float*,float*,float*,float*,float*);
__global__ void kernel_guassianintegrateblobs(int,int,int, float*,float*,float*, float*,float*,float*,float*,float*,float*);



#endif