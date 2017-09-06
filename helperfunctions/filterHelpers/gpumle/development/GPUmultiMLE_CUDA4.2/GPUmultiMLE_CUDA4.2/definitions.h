#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#define IMSZ 11
#define IMSZBIG 21
#define NK 128 //number of blocks to run in each kernel
#define pi 3.141592f
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#define max(a,b)            (((a) > (b)) ? (a) : (b))

// Thread block size GPUCRLB.cu
#define BSZ_CRLB 10
#define MEM_CRLB 60
#define MATMEM 1240
#define NUM_FIT_CRLB 5 //max number of multifit in previous step

// Thread block size GPUgenerateblobs.cu
#define BSZ_BLOBS 64
#define MEM_BLOBS 70
#define imMEM 4000

//GPUmultifit.cu
#define BSZ_FIT 64
#define MEM_FIT 3872
#define NUM_FIT 8

// Copied from GPUgaussMLEv2
#define BSZ 64			//!< max number of threads per block 
#define MEM 3872		//!< not used
//#define IMSZ 11			//!< not used
//#define IMSZBIG 21		//!< maximum fitting window size
//#define NK 128			//!< number of blocks to run in each kernel
//#define pi 3.141592f	//!< ensure a consistent value for pi
#define NV_P 4			//!< number of fitting parameters for MLEfit (x,y,bg,I)
#define NV_PS 5			//!< number of fitting parameters for MLEFit_sigma (x,y,bg,I,Sigma)
#define NV_PZ 5			//!< not used (x,y,bg,I,z)
#define NV_PS2 6		//!< number of fitting parameters for MLEFit_sigmaxy (x,y,bg,I,Sx,Sy)


#endif