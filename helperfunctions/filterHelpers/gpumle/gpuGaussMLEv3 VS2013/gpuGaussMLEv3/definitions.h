/*!
 * \file definitions.h
 * \author <your name>
 * \date April 20, 2012
 * \brief  If you are defining any constants they better be in this file and only this file.
 */

#ifndef DEFINITIONS_H
#define DEFINITIONS_H
#include <cuda_runtime.h>

#define BSZ 64			//!< max number of threads per block 
#define MEM 3872		//!< not used
#define IMSZ 11			//!< not used
#define IMSZBIG 21		//!< maximum fitting window size
#define NK 128			//!< number of blocks to run in each kernel
#define pi 3.141592f	//!< ensure a consistent value for pi
//#define BGMODEL 0 //!<backgroundmodel linear 1) gradient or  0) uniform.

//// MLERatio
#define NV_RH0 1		//!< number of fitting parameters for MLERatio (bg)
#define NV_RH1 2		//!< number of fitting parameters for MLERatio (bg,I)
#define Q_R 6		//!< number of quality parameters for MLERatio

//// MLEfit
#define NV_PH0 1			//!< number of fitting parameters for MLEfit (bg)
#define NV_PH1 4			//!< number of fitting parameters for MLEfit (x,y,bg,I)
#define Q_P 6

//MLEFit_sigma
#define NV_PSH0 1			//!< number of fitting parameters for MLEFit_sigma (bg)
#define NV_PSH1 5			//!< number of fitting parameters for MLEFit_sigma (x,y,bg,I,Sigma)
#define Q_PS 6


//// MLEfit3D
#define NV_PH03 1			//!< number of fitting parameters for MLEfit (bg)
#define NV_PH13 5			//!< number of fitting parameters for MLEfit (x,y,bg,I,z)
#define Q_P3 6

#define NV_PZH0 1			//!< not used (bg)
#define NV_PZH1 5			//!< not used (x,y,bg,I,z)
#define Q_PZ 6

//// MLEfitSigmaxy
#define NV_PS2H0 1 //!< number of fitting parameters for MLEfit (bg)
#define NV_PS2H1 6 //!< number of fitting parameters for MLEfit (x,y,bg,I,Sx,Sy)
#define Q_PS2 6

//// MLEfit3DSigma
#define NV_PH05 1			//!< number of fitting parameters for MLEfit (bg)
#define NV_PH15 7			//!< number of fitting parameters for MLEfit (x,y,z,bg,I,Sx,Sz)
#define Q_P5 6


//const numberOfVarFitType lookupNumberOfVarFitType[] = { { NV_RH0, NV_RH1, Q_R }, // (bg),(bg,I)
//{ NV_PH0, NV_PH1, Q_P }, // (bg),(x,y,bg,I)
//{ NV_PSH0, NV_PSH1, Q_PS }, // (bg),(x,y,bg,I,Sigma)
//{ NV_PH03, NV_PH13, Q_P3 }, // (bg), (x,y,bg,I,z)
//{ NV_PS2H0, NV_PS2H1, Q_PS2 }, // (bg),(x,y,bg,I,Sx,Sy)
//{ NV_PH05, NV_PH15, Q_P5 }, //  (bg)(x,y,z,bg,I,Sx,Sz)
//};

//


//
////MLEFit_sigma
//#define NV_PSH0 1			//!< number of fitting parameters for MLEFit_sigma (bg)
//#define NV_PSH1 5			//!< number of fitting parameters for MLEFit_sigma (x,y,bg,I,Sigma)
//#define Q_PS 6
//
////not used
//#define NV_PZH0 1			//!< not used (bg)
//#define NV_PZH1 5			//!< not used (x,y,bg,I,z)
//#define Q_PZ 6
//
////MLEFit_sigmaxy
//#define NV_PS2H0 1		//!< number of fitting parameters for MLEFit_sigmaxy (bg)
//#define NV_PS2H1 6		//!< number of fitting parameters for MLEFit_sigmaxy (x,y,bg,I,Sx,Sy)
//#define Q_PS2 6

typedef struct
{
	//Device buffers
	float *d_Data;
	float *h_Data;

	int BlockSize;
	int NKS;

	size_t maxNumberOfFits;

	size_t MEMS;
	
	size_t MAXMEMS;
	
	int IMBIGS;

	int Nfits;
	float gpuDist;

	float* d_Parameters;
	int szParameters;
	float* d_CRLBs;
	int szCRLBs;
	float* d_LogLikelihood;
	int szLogLikelihood;
	int maxGrid;

	float* h_Parameters;
	float* h_CRLBs;
	float* h_LogLikelihood;

	//Stream for asynchronous command execution
	cudaStream_t stream;
	int dimBlock;
	int dimGrid;
	int major;
	char name[256];

} TGPUplan;

//void CUDAERROR(const char *instr,int lineNumber);
//void cudasafe( cudaError_t err, char* str, int lineNumber);

#endif