/*!
 * \file MatPropLib.cuh
 * \author Keith Lidke & Carlas Smith
 * \date May 14, 2014
 * \brief This provides the probability helper functions.
 */
#include "PropLib.h"
#define _USE_MATH_DEFINES
#include "math.h"


__device__ float normcdf(float x, float mu, float sigma){
     return 0.5*(1+erf((x-mu)/(sqrt(2.0)*sigma)));   
}

__device__ float norminv(float p, float mu, float sigma){
     return  mu+sigma*sqrt(2.0)*erfinv(2.0*p-1);
}
__device__ float normpdf(float x, float mu, float sigma){
     return  1/(sigma*sqrt(2.0*M_PI))*exp(-pow(x-mu,2)/(2.0*pow(sigma,2)));
}

    