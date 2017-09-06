#ifndef IMAGE_OPERATION_H
#define IMAGE_OPERATION_H

void MakeSubregions(const int candicc, const float *data, const int szinX, const int szinY, const int szinZ, float *candx,
					float *candy, float *candz, float boxsz, float *dataout, float *left, float *top);
void GPUmultifit(int Nfitraw, float *in_sub,const int sz, float psfsigma,const int iterations,int fitnum, float Nave, 
				 const float Nsigmafit,float threshold,float *xx,float *yy,float *nn,float *bb,float *div);
void GPUCRLB(int candicc,const int boxsz,float psfsigma,int fitnum, float *x, float *y, float *n, 
			 float *bg, float Nave,float Nsigmafit,float *CRLBarray, float *covariance);
void GPUgenerateblobs(int blobn,int xsz,int ysz,float *xx,float *yy, float *nn, float *LAx, float *LAy, float *cov, float *im);

#endif