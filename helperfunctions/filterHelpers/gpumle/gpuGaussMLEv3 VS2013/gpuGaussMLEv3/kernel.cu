#include <cuda_runtime.h>
#include "definitions.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MatrixLib.cuh"
#include "GPUgaussLib.cuh"

//*******************************************************************************************
//theta is: {N,bg}
__global__ void kernel_MLERatio(const float *d_data, float PSFSigma, int sz, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits){
	/*!
	* \brief basic MLE fitting kernel.  No additional parameters are computed.
	* \param d_data array of subregions to fit copied to GPU
	* \param PSFSigma sigma of the point spread function
	* \param sz nxn size of the subregion to fit
	* \param iterations number of iterations for solution to converge
	* \param d_Parameters array of fitting parameters to return for each subregion
	* \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	* \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	* \param Nfits number of subregions to fit
	*/
	float M[NV_RH1*NV_RH1], Diag[NV_RH1], Minv[NV_RH1*NV_RH1];
	int tx = threadIdx.x;
	int bx = blockIdx.x;
	int BlockSize = blockDim.x;
	int ii, jj, kk, ll;
	float model, cf, df, data;
	float Div;
	float PSFy, PSFx;
	float dudt[NV_RH1];
	float d2udt2[NV_RH1];
	float NR_Numerator[NV_RH1], NR_Denominator[NV_RH1];
	float thetaH1[NV_RH1];
	float thetaH0[NV_RH0];
	//float maxjump[NV_RH1] = { 1e2f, 2e0f };
	//float gamma[NV_RH1] = { 0.5f, 1.0f };
	float Nmax;
	float logModel;
	float LR[Q_R];

	//Prevent read/write past end of array
	if ((bx*BlockSize + tx) >= Nfits) return;

	memset(M, 0, NV_RH1*NV_RH1*sizeof(float));
	memset(Minv, 0, NV_RH1*NV_RH1*sizeof(float));

	//load data
	const float *s_data = d_data + (sz*sz*bx*BlockSize + sz*sz*tx);

	kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &thetaH1[1]);
	thetaH1[0] = max(0.1f, (Nmax - thetaH1[1]) * 4 * pi*PSFSigma*PSFSigma);

	for (kk = 0; kk<iterations; kk++) {//main iterative loop

		//initialize
		memset(NR_Numerator, 0, NV_RH1*sizeof(float));
		memset(NR_Denominator, 0, NV_RH1*sizeof(float));

		for (ii = 0; ii<sz; ii++) for (jj = 0; jj<sz; jj++) {
			PSFx = kernel_IntGauss1D(ii, (sz-1) / 2.0, PSFSigma);
			PSFy = kernel_IntGauss1D(jj, (sz-1) / 2.0, PSFSigma);

			model = thetaH1[1] + thetaH1[0] * PSFx*PSFy;
			data = s_data[sz*jj + ii];

			//calculating derivatives
			dudt[0] = PSFx*PSFy;
			d2udt2[0] = 0.0f;
			dudt[1] = 1.0f;
			d2udt2[1] = 0.0f;

			cf = 0.0f;
			df = 0.0f;
			if (model>10e-3f) cf = data / model - 1;
			if (model>10e-3f) df = data / pow(model, 2);
			cf = min(cf, 10e4f);
			df = min(df, 10e4f);

			for (ll = 0; ll<NV_RH1; ll++){
				NR_Numerator[ll] += dudt[ll] * cf;
				NR_Denominator[ll] += d2udt2[ll] * cf - pow(dudt[ll], 2)*df;
			}
		}

		// Any other constraints
		thetaH1[0] -= min(max(NR_Numerator[0] / NR_Denominator[0] / 2.0, -thetaH1[0]), thetaH1[0] / 2.0);
		thetaH1[0] = max(thetaH1[0], Nmax/2.0f);

		thetaH1[1] -= NR_Numerator[1] / NR_Denominator[1];
		thetaH1[1] = max(thetaH1[1], 0.01f);

	}

	// ML estimate of background model
	thetaH0[0] = 0.0;
	for (ii = 0; ii<sz; ii++) for (jj = 0; jj<sz; jj++) {
		thetaH0[0] += s_data[sz*jj + ii];
	}
	thetaH0[0] = thetaH0[0] / pow((float)sz, 2);

	// Calculating the CRLB and LogLikelihoodRatio
	Div = 0.0;
	for (ii = 0; ii<sz; ii++) for (jj = 0; jj<sz; jj++) {
		PSFx = kernel_IntGauss1D(ii, (sz - 1) / 2.0, PSFSigma);
		PSFy = kernel_IntGauss1D(jj, (sz - 1) / 2.0, PSFSigma);

		model = thetaH1[1] + thetaH1[0] * PSFx*PSFy;
		data = s_data[sz*jj + ii];

		//calculating derivatives
		dudt[0] = PSFx*PSFy;
		dudt[1] = 1.0f;

		//Building the Fisher Information Matrix
		for (kk = 0; kk<NV_RH1; kk++)for (ll = kk; ll<NV_RH1; ll++){
			M[kk*NV_RH1 + ll] += dudt[ll] * dudt[kk] / model;
			M[ll*NV_RH1 + kk] = M[kk*NV_RH1 + ll];
		}

		//LogLikelyhood
		logModel = model / (thetaH0[0] + 1e-5);
		if (logModel>0 && data > 0)
			Div += 2 * (data*log(logModel + 1e-5) - model + thetaH0[0]);
	}

	// Matrix inverse (CRLB=F^-1) and output assigments
	kernel_MatInvN(M, Minv, Diag, NV_RH1);
	kernel_CalcLLRProp(Diag[0], thetaH1[0], Div, LR);

	//write to global arrays
	for (kk = 0; kk<NV_RH1; kk++)
		d_Parameters[kk + (NV_RH1 + NV_RH0)*(BlockSize*bx + tx)] = thetaH1[kk];
	for (kk = 0; kk<NV_RH0; kk++)
		d_Parameters[(NV_RH1 + kk) + (NV_RH1 + NV_RH0)*(BlockSize*bx + tx)] = thetaH0[kk];
	for (kk = 0; kk<Q_R; kk++)
		d_LogLikelihood[kk + Q_R * (BlockSize*bx + tx)] = LR[kk];
	for (kk = 0; kk<NV_RH1; kk++)
		d_CRLBs[kk + NV_RH1*(BlockSize*bx + tx)] = Diag[kk];
	return;
}
//*******************************************************************************************
//theta is: {x,y,N,bg}
__global__ void kernel_MLEFit(const float *d_data, float PSFSigma, int sz, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits){
	/*!
	* \brief basic MLE fitting kernel.  No additional parameters are computed.
	* \param d_data array of subregions to fit copied to GPU
	* \param PSFSigma sigma of the point spread function
	* \param sz nxn size of the subregion to fit
	* \param iterations number of iterations for solution to converge
	* \param d_Parameters array of fitting parameters to return for each subregion
	* \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	* \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	* \param Nfits number of subregions to fit
	*/
	float M[NV_PH1*NV_PH1], Diag[NV_PH1], Minv[NV_PH1*NV_PH1];
	int tx = threadIdx.x;
	int bx = blockIdx.x;
	int BlockSize = blockDim.x;
	int ii, jj, kk, ll;
	float model, cf, df, data;
	float Div;
	float PSFy, PSFx;
	//int NV = NV_PH1;
	float dudt[NV_PH1];
	float d2udt2[NV_PH1];
	float NR_Numerator[NV_PH1], NR_Denominator[NV_PH1];
	float thetaH1[NV_PH1];
	float thetaH0[NV_PH0];
	float maxjump[NV_PH1] = { 1e0f, 1e0f, 1e2f, 2e0f };
	float gamma[NV_PH1] = { 1.0f, 1.0f, 0.5f, 1.0f };
	float LR[Q_P];
	float Nmax;
	float logModel;

	//Prevent read/write past end of array
	if ((bx*BlockSize + tx) >= Nfits) return;

	memset(M, 0, NV_PH1*NV_PH1*sizeof(float));
	memset(Minv, 0, NV_PH1*NV_PH1*sizeof(float));

	const float *s_data = d_data + (sz*sz*bx*BlockSize + sz*sz*tx);

	//initial values
	kernel_CenterofMass2D(sz, s_data, &thetaH1[0], &thetaH1[1], 0, 0);
	kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &thetaH1[3]);
	thetaH1[2] = max(0.1f, (Nmax - thetaH1[3]) * 2 * pi*PSFSigma*PSFSigma);

	for (kk = 0; kk<iterations; kk++) {//main iterative loop

		memset(NR_Numerator, 0, NV_PH1*sizeof(float));
		memset(NR_Denominator, 0, NV_PH1*sizeof(float));

		for (ii = 0; ii<sz; ii++) for (jj = 0; jj<sz; jj++) {
			PSFx = kernel_IntGauss1D(ii, thetaH1[0], PSFSigma);
			PSFy = kernel_IntGauss1D(jj, thetaH1[1], PSFSigma);

			model = thetaH1[3] + thetaH1[2] * PSFx*PSFy;
			data = s_data[sz*jj + ii];

			//calculating derivatives
			kernel_DerivativeIntGauss1D(ii, thetaH1[0], PSFSigma, thetaH1[2], PSFy, &dudt[0], &d2udt2[0]);
			kernel_DerivativeIntGauss1D(jj, thetaH1[1], PSFSigma, thetaH1[2], PSFx, &dudt[1], &d2udt2[1]);
			dudt[2] = PSFx*PSFy;
			d2udt2[2] = 0.0f;
			dudt[3] = 1.0f;
			d2udt2[3] = 0.0f;

			cf = 0.0f;
			df = 0.0f;
			if (model>10e-3f) cf = data / model - 1;
			if (model>10e-3f) df = data / pow(model, 2);
			cf = min(cf, 10e4f);
			df = min(df, 10e4f);

			for (ll = 0; ll<NV_PH1; ll++){
				NR_Numerator[ll] += dudt[ll] * cf;
				NR_Denominator[ll] += d2udt2[ll] * cf - pow(dudt[ll], 2)*df;
			}
		}

		// The update
		if (kk<2)
			for (ll = 0; ll<NV_PH1; ll++)
				thetaH1[ll] -= gamma[ll] * min(max(NR_Numerator[ll] / NR_Denominator[ll], -maxjump[ll]), maxjump[ll]);
		else
			for (ll = 0; ll<NV_PH1; ll++)
				thetaH1[ll] -= min(max(NR_Numerator[ll] / NR_Denominator[ll], -maxjump[ll]), maxjump[ll]);

		// Any other constraints
		thetaH1[2] = max(thetaH1[2], 1.0f);
		thetaH1[3] = max(thetaH1[3], 0.01f);

	}

	//Estimate background model   
	thetaH0[0] = 0.0;
	for (ii = 0; ii<sz; ii++) for (jj = 0; jj<sz; jj++) {
		thetaH0[0] += s_data[sz*jj + ii];
	}
	thetaH0[0] = thetaH0[0] / pow((float)sz, 2);
	// Calculating the CRLB and LogLikelihood
	Div = 0.0;
	for (ii = 0; ii<sz; ii++) for (jj = 0; jj<sz; jj++) {
		PSFx = kernel_IntGauss1D(ii, thetaH1[0], PSFSigma);
		PSFy = kernel_IntGauss1D(jj, thetaH1[1], PSFSigma);

		model = thetaH1[3] + max(thetaH1[2], thetaH1[3])*PSFx*PSFy;
		data = s_data[sz*jj + ii];

		//calculating derivatives
		kernel_DerivativeIntGauss1D(ii, thetaH1[0], PSFSigma, thetaH1[2], PSFy, &dudt[0], NULL);
		kernel_DerivativeIntGauss1D(jj, thetaH1[1], PSFSigma, thetaH1[2], PSFx, &dudt[1], NULL);
		dudt[2] = PSFx*PSFy;
		dudt[3] = 1.0f;

		//Building the Fisher Information Matrix
		for (kk = 0; kk<NV_PH1; kk++)for (ll = kk; ll<NV_PH1; ll++){
			M[kk*NV_PH1 + ll] += dudt[ll] * dudt[kk] / model;
			M[ll*NV_PH1 + kk] = M[kk*NV_PH1 + ll];
		}

		//LogLikelyhood
		logModel = model / (thetaH0[0] + 1e-5);
		if (logModel>0 && data > 0)
			Div += 2 * (data*log(logModel + 1e-5) - model + thetaH0[0]);
	}

	// Matrix inverse (CRLB=F^-1) and output assigments
	kernel_MatInvN(M, Minv, Diag, NV_PH1);
	kernel_CalcLLRProp(Diag[2], thetaH1[2], Div, LR);

	for (kk = 0; kk<NV_PH1; kk++) 
		d_Parameters[kk + (NV_PH1 + NV_PH0)*(BlockSize*bx + tx)] = thetaH1[kk];
	for (kk = 0; kk<NV_PH0; kk++)
		d_Parameters[(NV_PH1 + kk) + (NV_PH1 + NV_PH0)*(BlockSize*bx + tx)] = thetaH0[kk];
	for (kk = 0; kk<Q_P; kk++) 
		d_LogLikelihood[kk + Q_P * (BlockSize*bx + tx)] = LR[kk];
	for (kk = 0; kk<NV_PH1; kk++) 
		d_CRLBs[kk + NV_PH1*(BlockSize*bx + tx)] = Diag[kk];

	return;
}

//*******************************************************************************************
__global__ void kernel_MLEFitSigma(float *d_data, float PSFSigma, int sz, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits){
		/*!
		* \brief basic MLE fitting kernel.  No additional parameters are computed.
		* \param d_data array of subregions to fit copied to GPU
		* \param PSFSigma sigma of the point spread function
		* \param sz nxn size of the subregion to fit
		* \param iterations number of iterations for solution to converge
		* \param d_Parameters array of fitting parameters to return for each subregion
		* \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
		* \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
		* \param Nfits number of subregions to fit
		*/

		float M[NV_PSH1*NV_PSH1], Diag[NV_PSH1], Minv[NV_PSH1*NV_PSH1];
		int tx = threadIdx.x;
		int bx = blockIdx.x;
		int BlockSize = blockDim.x;
		int ii, jj, kk, ll;
		float model, cf, df, data;
		float Div;
		float PSFy, PSFx;
		float dudt[NV_PSH1];
		float d2udt2[NV_PSH1];
		float NR_Numerator[NV_PSH1], NR_Denominator[NV_PSH1];
		float thetaH1[NV_PSH1];
		float thetaH0[NV_PSH0];

		float maxjump[NV_PSH1] = { 1e0f, 1e0f, 1e2f, 2e0f, 5e-1f };
		float gamma[NV_PSH1] = { 1.0f, 1.0f, 0.5f, 1.0f, 1.0f };
		float Nmax;
		float logModel;
		float LR[Q_PS];

		//Prevent read/write past end of array
		if ((bx*BlockSize + tx) >= Nfits) return;

		memset(M, 0, NV_PSH1*NV_PSH1*sizeof(float));
		memset(Minv, 0, NV_PSH1*NV_PSH1*sizeof(float));
		//load data
		const float *s_data = d_data + (sz*sz*bx*BlockSize + sz*sz*tx);

		//initial values
		kernel_CenterofMass2D(sz, s_data, &thetaH1[0], &thetaH1[1], 0, 0);
		kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &thetaH1[3]);
		thetaH1[2] = max(0.1f, (Nmax - thetaH1[3]) * 2 * pi*PSFSigma*PSFSigma);
		thetaH1[4] = PSFSigma;


		for (kk = 0; kk<iterations; kk++) {//main iterative loop

			//initialize
			memset(NR_Numerator, 0, NV_PSH1*sizeof(float));
			memset(NR_Denominator, 0, NV_PSH1*sizeof(float));

			for (ii = 0; ii<sz; ii++) for (jj = 0; jj<sz; jj++) {
				PSFx = kernel_IntGauss1D(ii, thetaH1[0], thetaH1[4]);
				PSFy = kernel_IntGauss1D(jj, thetaH1[1], thetaH1[4]);

				model = thetaH1[3] + thetaH1[2] * PSFx*PSFy;
				data = s_data[sz*jj + ii];

				//calculating derivatives
				kernel_DerivativeIntGaussPSF1D(ii, thetaH1[0], thetaH1[4], thetaH1[2], PSFy, &dudt[0], &d2udt2[0]);
				kernel_DerivativeIntGaussPSF1D(jj, thetaH1[1], thetaH1[4], thetaH1[2], PSFx, &dudt[1], &d2udt2[1]);
				kernel_DerivativeIntGaussPSF2DSigma(ii, jj, thetaH1[0], thetaH1[1], thetaH1[4], thetaH1[2], PSFx, PSFy, &dudt[4], &d2udt2[4]);
				dudt[2] = PSFx*PSFy;
				d2udt2[2] = 0.0f;
				dudt[3] = 1.0f;
				d2udt2[3] = 0.0f;

				cf = 0.0f;
				df = 0.0f;
				if (model>10e-3f) cf = data / model - 1;
				if (model>10e-3f) df = data / pow(model, 2);
				cf = min(cf, 10e4f);
				df = min(df, 10e4f);

				for (ll = 0; ll<NV_PSH1; ll++){
					NR_Numerator[ll] += dudt[ll] * cf;
					NR_Denominator[ll] += d2udt2[ll] * cf - pow(dudt[ll], 2)*df;
				}
			}

			// The update
			if (kk<5)
				for (ll = 0; ll<NV_PSH1; ll++)
					thetaH1[ll] -= gamma[ll] * min(max(NR_Numerator[ll] / NR_Denominator[ll], -maxjump[ll]), maxjump[ll]);
			else
				for (ll = 0; ll<NV_PSH1; ll++)
					thetaH1[ll] -= min(max(NR_Numerator[ll] / NR_Denominator[ll], -maxjump[ll]), maxjump[ll]);

			// Any other constraints
			thetaH1[2] = max(thetaH1[2], 1.0f);
			thetaH1[3] = max(thetaH1[3], 0.01f);
			thetaH1[4] = max(thetaH1[4], 0.5f);
			thetaH1[4] = min(thetaH1[4], sz / 2.0f);
		}

		//Estimate background model   
		thetaH0[0] = 0.0;
		for (ii = 0; ii<sz; ii++) for (jj = 0; jj<sz; jj++) {
			thetaH0[0] += s_data[sz*jj + ii];
		}
		thetaH0[0] = thetaH0[0] / pow((float)sz, 2);

		// Calculating the CRLB and LogLikelihood
		Div = 0.0f;
		for (ii = 0; ii<sz; ii++) for (jj = 0; jj<sz; jj++) {
			PSFx = kernel_IntGauss1D(ii, thetaH1[0], PSFSigma);
			PSFy = kernel_IntGauss1D(jj, thetaH1[1], PSFSigma);

			model = thetaH1[3] + thetaH1[2] * PSFx*PSFy;
			data = s_data[sz*jj + ii];

			//calculating derivatives
			kernel_DerivativeIntGaussPSF1D(ii, thetaH1[0], thetaH1[4], thetaH1[2], PSFy, &dudt[0], NULL);
			kernel_DerivativeIntGaussPSF1D(jj, thetaH1[1], thetaH1[4], thetaH1[2], PSFx, &dudt[1], NULL);
			kernel_DerivativeIntGaussPSF2DSigma(ii, jj, thetaH1[0], thetaH1[1], thetaH1[4], thetaH1[2], PSFx, PSFy, &dudt[4], NULL);
			dudt[2] = PSFx*PSFy;
			dudt[3] = 1.0f;

			//Building the Fisher Information Matrix
			for (kk = 0; kk<NV_PSH1; kk++)for (ll = kk; ll<NV_PSH1; ll++){
				M[kk*NV_PSH1 + ll] += dudt[ll] * dudt[kk] / model;
				M[ll*NV_PSH1 + kk] = M[kk*NV_PSH1 + ll];
			}

			//LogLikelyhood
			logModel = model / (thetaH0[0] + 1e-5);
			if (logModel>0 && data > 0)
				Div += 2 * (data*log(logModel + 1e-5) - model + thetaH0[0]);
		}

		// Matrix inverse (CRLB=F^-1) and output assigments
		kernel_MatInvN(M, Minv, Diag, NV_PSH1);
		kernel_CalcLLRProp(Diag[2], thetaH1[2], Div, LR);

		for (kk = 0; kk<NV_PSH1; kk++)
			d_Parameters[kk + (NV_PSH1 + NV_PSH0)*(BlockSize*bx + tx)] = thetaH1[kk];
		for (kk = 0; kk<NV_PSH0; kk++)
			d_Parameters[(NV_PSH1 + kk) + (NV_PSH1 + NV_PSH0)*(BlockSize*bx + tx)] = thetaH0[kk];
		for (kk = 0; kk<Q_PS; kk++)
			d_LogLikelihood[kk + Q_PS * (BlockSize*bx + tx)] = LR[kk];
		for (kk = 0; kk<NV_PSH1; kk++)
			d_CRLBs[kk + NV_PSH1*(BlockSize*bx + tx)] = Diag[kk];

		return;
	}



__global__ void kernel_MLEFit3D(const float *d_data, float PSFSigmax, float PSFSigmaz, int szXY, int szZ, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits){
	/*!
	* \brief basic MLE fitting kernel.  No additional parameters are computed.
	* \param d_data array of subregions to fit copied to GPU
	* \param PSFSigma sigma of the point spread function
	* \param sz nxn size of the subregion to fit
	* \param iterations number of iterations for solution to converge
	* \param d_Parameters array of fitting parameters to return for each subregion
	* \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	* \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	* \param Nfits number of subregions to fit
	*/
	
	//Dynamic allocation is slow. Maybe this will change overtime but right now keep it static!
	float M[NV_PH13*NV_PH13], Diag[NV_PH13], Minv[NV_PH13*NV_PH13];
	int tx = threadIdx.x;
	int bx = blockIdx.x;
	int BlockSize = blockDim.x;
	int hh, ii, jj, kk, ll;
	float model, cf, df, data;
	float Div;
	float PSFy, PSFx, PSFz;
	float *dudt = new float[NV_PH13 + szZ - 1];
	float *d2udt2 = new float[NV_PH13 + szZ - 1];

	float NR_Numerator[NV_PH13], NR_Denominator[NV_PH13];
	float *thetaH1 = new float[NV_PH13 + szZ - 1];	
	float *thetaH0 = new float[NV_PH03+szZ-1];

	// last value is used for all variables bigger then NV_PH13
	float maxjump[NV_PH13] = { 1e0f, 1e0f, 1e0f, 1e2f, 2e0f }; 
	float gamma[NV_PH13] = { 1.0f, 1.0f, 1.0f, 0.5f, 1.0f };
	float LR[Q_P3];
	float Nmax;
	float logModel;

	//Prevent read/write past end of array
	if ((bx*BlockSize + tx) >= Nfits) return;

	memset(M, 0, NV_PH13*NV_PH13*sizeof(float));
	memset(Minv, 0, NV_PH13*NV_PH13*sizeof(float));

	const float *s_data = d_data + (szXY*szXY*szZ*bx*BlockSize + szXY*szXY*szZ*tx);

	////initial values
	kernel_CenterofMass3D(szXY, szZ, s_data, &thetaH1[0], &thetaH1[1], &thetaH1[2], 0, 0, 0);
	kernel_GaussFMaxMin3D(szXY, szZ, PSFSigmax, PSFSigmaz, s_data, &Nmax, &thetaH1[4]);

	thetaH1[3] = max(0.1f, (Nmax - thetaH1[4]) * pow(sqrt(2 * pi),3) * PSFSigmax * PSFSigmax * PSFSigmaz);
	for (kk = 0; kk<iterations; kk++) {//main iterative loop
		
		memset(NR_Numerator, 0, NV_PH13*sizeof(float));
		memset(NR_Denominator, 0, NV_PH13*sizeof(float));

		for (hh = 0; hh < szZ; hh++) for (ii = 0; ii < szXY; ii++) for (jj = 0; jj < szXY; jj++) {
			PSFx = kernel_IntGauss1D(ii, thetaH1[0], PSFSigmax);
			PSFy = kernel_IntGauss1D(jj, thetaH1[1], PSFSigmax);
			PSFz = kernel_IntGauss1D(hh, thetaH1[2], PSFSigmaz);

			model = thetaH1[4+hh] + thetaH1[3] * PSFx*PSFy*PSFz;
			data = s_data[hh*szXY*szXY +szXY*ii+jj];

			//calculating derivatives
			kernel_DerivativeIntGaussPSF1D(ii, thetaH1[0], PSFSigmax, thetaH1[3], PSFy, PSFz, &dudt[0], &d2udt2[0]);
			kernel_DerivativeIntGaussPSF1D(jj, thetaH1[1], PSFSigmax, thetaH1[3], PSFx, PSFz, &dudt[1], &d2udt2[1]);
			kernel_DerivativeIntGaussPSF1D(hh, thetaH1[2], PSFSigmaz, thetaH1[3], PSFx, PSFy, &dudt[2], &d2udt2[2]);
			dudt[3] = PSFx*PSFy*PSFz;
			d2udt2[3] = 0.0f;
			dudt[4+hh] = 1.0f;
			d2udt2[4+hh] = 0.0f;

			cf = 0.0f;
			df = 0.0f;
			if (model>10e-3f) cf = data / model - 1;
			if (model>10e-3f) df = data / pow(model, 2);
			cf = min(cf, 10e4f);
			df = min(df, 10e4f);

			for (ll = 0; ll<NV_PH13; ll++){
				NR_Numerator[ll] += dudt[ll] * cf;
				NR_Denominator[ll] += d2udt2[ll] * cf - pow(dudt[ll], 2)*df;
			}
		}

		//// The update
		if (kk<3)
			for (ll = 0; ll<NV_PH13 + szZ - 1; ll++)
				thetaH1[ll] -= gamma[min(NV_PH13, ll)] * min(max(NR_Numerator[ll] / NR_Denominator[ll], -maxjump[min(NV_PH13, ll)]), maxjump[min(NV_PH13, ll)]);
		else
			for (ll = 0; ll<NV_PH13 + szZ - 1; ll++)
				thetaH1[ll] -= min(max(NR_Numerator[ll] / NR_Denominator[ll], -maxjump[min(NV_PH13, ll)]), maxjump[min(NV_PH13, ll)]);

		// Any other constraints
		thetaH1[3] = max(thetaH1[3], 1.0f);
		for (hh = 4; hh<NV_PH13 + szZ - 1; hh++)
			thetaH1[hh] = max(thetaH1[hh], 0.1f);
	}

	//Estimate background model   
	thetaH0[0] = 0.0;
	for (hh = 0; hh < szZ; hh++){
		for (ii = 0; ii < szXY; ii++) for (jj = 0; jj < szXY; jj++) {
			//TODO add oo estimation!!
			thetaH0[hh] += s_data[hh*szXY*szXY + szXY*ii + jj];
		}
		thetaH0[hh] = thetaH0[hh] / pow((float)szXY, 2);
	}
	// Calculating the CRLB and LogLikelihood
	Div = 0.0;
	for (hh = 0; hh < szZ; hh++) for (ii = 0; ii<szXY; ii++) for (jj = 0; jj<szXY; jj++) {
		PSFx = kernel_IntGauss1D(ii, thetaH1[0], PSFSigmax);
		PSFy = kernel_IntGauss1D(jj, thetaH1[1], PSFSigmax);
		PSFz = kernel_IntGauss1D(hh, thetaH1[2], PSFSigmaz);

		model = thetaH1[4+hh] + max(thetaH1[3], thetaH1[4+hh])*PSFx*PSFy*PSFz;
		data = s_data[hh*szXY*szXY +szXY*ii+jj];

		//calculating derivatives
		kernel_DerivativeIntGaussPSF1D(ii, thetaH1[0], PSFSigmax, thetaH1[3], PSFy, PSFz, &dudt[0], NULL);
		kernel_DerivativeIntGaussPSF1D(jj, thetaH1[1], PSFSigmax, thetaH1[3], PSFx, PSFz, &dudt[1], NULL);
		kernel_DerivativeIntGaussPSF1D(hh, thetaH1[2], PSFSigmaz, thetaH1[3], PSFx, PSFy, &dudt[2], NULL);
		dudt[3] = PSFx*PSFy*PSFz;
		dudt[4] = 1.0f;

		//Building the Fisher Information Matrix
		for (kk = 0; kk<(NV_PH13 + szZ - 1); kk++)for (ll = kk; ll<(NV_PH13 + szZ - 1); ll++){
			M[kk*(NV_PH13 + szZ - 1) + ll] += dudt[ll] * dudt[kk] / model;
			M[ll*(NV_PH13 + szZ - 1) + kk] = M[kk*(NV_PH13 + szZ - 1) + ll];
		}

		//LogLikelyhood
		logModel = model / (thetaH0[hh] + 1e-5);
		if (logModel>0 && data > 0)
			Div += 2 * (data*log(logModel + 1e-5) - model + thetaH0[hh]);
	}

	//// Matrix inverse (CRLB=F^-1) and output assigments
	kernel_MatInvN(M, Minv, Diag, (NV_PH13 + szZ - 1));
	kernel_CalcLLRProp(Diag[3], thetaH1[3], Div, LR);

	for (kk = 0; kk<(NV_PH13 + szZ - 1); kk++)
		d_Parameters[kk + ((NV_PH13 + szZ - 1) + (NV_PH03 + szZ - 1))*(BlockSize*bx + tx)] = thetaH1[kk];
	for (kk = 0; kk<(NV_PH03 + szZ - 1); kk++)
		d_Parameters[((NV_PH13 + szZ - 1) + kk) + ((NV_PH13 + szZ - 1) + (NV_PH03 + szZ - 1))*(BlockSize*bx + tx)] = thetaH0[kk];
	for (kk = 0; kk<Q_P3; kk++)
		d_LogLikelihood[kk + Q_P3 * (BlockSize*bx + tx)] = LR[kk];
	for (kk = 0; kk<(NV_PH13 + szZ - 1); kk++)
		d_CRLBs[kk + (NV_PH13 + szZ - 1)*(BlockSize*bx + tx)] = Diag[kk];

	return;
}

__global__ void kernel_MLEFit3DSigma(const float *d_data, float PSFSigmax, float PSFSigmaz, int szXY, int szZ, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits){
	/*!
	* \brief basic MLE fitting kernel.  No additional parameters are computed.
	* \param d_data array of subregions to fit copied to GPU
	* \param PSFSigma sigma of the point spread function
	* \param sz nxn size of the subregion to fit
	* \param iterations number of iterations for solution to converge
	* \param d_Parameters array of fitting parameters to return for each subregion
	* \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	* \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	* \param Nfits number of subregions to fit
	*/

	//Dynamic allocation is slow. Maybe this will change overtime but right now keep it static!
	float M[NV_PH15*NV_PH15], Diag[NV_PH15], Minv[NV_PH15*NV_PH15];
	int tx = threadIdx.x;
	int bx = blockIdx.x;
	int BlockSize = blockDim.x;
	int hh, ii, jj, kk, ll;
	float model, cf, df, data;
	float Div;
	float PSFy, PSFx, PSFz;
	float dudt[NV_PH15];
	float d2udt2[NV_PH15];
	float NR_Numerator[NV_PH15], NR_Denominator[NV_PH15];
	float thetaH1[NV_PH15];
	float thetaH0[NV_PH05];
	float maxjump[NV_PH15] = { 1e0f, 1e0f, 1e0f, 1e2f, 2e0f, 5e-1f, 5e-1f };
	float gamma[NV_PH15] = { 1.0f, 1.0f, 1.0f, 0.5f, 1.0f, 1.0f, 1.0f };
	float LR[Q_P5];
	float Nmax;
	float logModel;

	//Prevent read/write past end of array
	if ((bx*BlockSize + tx) >= Nfits) return;

	memset(M, 0, NV_PH15*NV_PH15*sizeof(float));
	memset(Minv, 0, NV_PH15*NV_PH15*sizeof(float));

	const float *s_data = d_data + (szXY*szXY*szZ*bx*BlockSize + szXY*szXY*szZ*tx);

	////initial values
	kernel_CenterofMass3D(szXY, szZ, s_data, &thetaH1[0], &thetaH1[1], &thetaH1[2], 0, 0, 0);
	kernel_GaussFMaxMin3D(szXY, szZ, PSFSigmax, PSFSigmaz, s_data, &Nmax, &thetaH1[4]);
	thetaH1[5] = PSFSigmax;
	thetaH1[6] = PSFSigmaz;

	thetaH1[3] = max(0.1f, (Nmax - thetaH1[4]) * pow(sqrt(2 * pi), 3) * PSFSigmax * PSFSigmax * PSFSigmaz);
	for (kk = 0; kk<iterations; kk++) {//main iterative loop

		memset(NR_Numerator, 0, NV_PH15*sizeof(float));
		memset(NR_Denominator, 0, NV_PH15*sizeof(float));

		for (hh = 0; hh < szZ; hh++) for (ii = 0; ii < szXY; ii++) for (jj = 0; jj < szXY; jj++) {
			PSFx = kernel_IntGauss1D(ii, thetaH1[0], thetaH1[5]);
			PSFy = kernel_IntGauss1D(jj, thetaH1[1], thetaH1[5]);
			PSFz = kernel_IntGauss1D(hh, thetaH1[2], thetaH1[6]);

			model = thetaH1[4] + thetaH1[3] * PSFx*PSFy*PSFz;
			data = s_data[hh*szXY*szXY + szXY*ii + jj];

			//calculating derivatives
			kernel_DerivativeIntGaussPSF1D(ii, thetaH1[0], thetaH1[5], thetaH1[3], PSFy, PSFz, &dudt[0], &d2udt2[0]);
			kernel_DerivativeIntGaussPSF1D(jj, thetaH1[1], thetaH1[5], thetaH1[3], PSFx, PSFz, &dudt[1], &d2udt2[1]);
			kernel_DerivativeIntGaussPSF1D(hh, thetaH1[2], thetaH1[6], thetaH1[3], PSFx, PSFy, &dudt[2], &d2udt2[2]);
			
			dudt[3] = PSFx*PSFy*PSFz;
			d2udt2[3] = 0.0f;
			dudt[4] = 1.0f;
			d2udt2[4] = 0.0f;

			kernel_DerivativeIntGaussPSF3DSigma(ii, jj, hh, thetaH1[0], thetaH1[1], thetaH1[2],
				thetaH1[5], thetaH1[6], thetaH1[3], PSFx, PSFy, PSFz, &dudt[5], &dudt[6], &d2udt2[5], &d2udt2[6]);

			cf = 0.0f;
			df = 0.0f;
			if (model>10e-3f) cf = data / model - 1;
			if (model>10e-3f) df = data / pow(model, 2);
			cf = min(cf, 10e4f);
			df = min(df, 10e4f);

			for (ll = 0; ll<NV_PH15; ll++){
				NR_Numerator[ll] += dudt[ll] * cf;
				NR_Denominator[ll] += d2udt2[ll] * cf - pow(dudt[ll], 2)*df;
			}
		}

		// The update
		if (kk<3)
			for (ll = 0; ll<NV_PH15; ll++)
				thetaH1[ll] -= gamma[ll] * min(max(NR_Numerator[ll] / NR_Denominator[ll], -maxjump[ll]), maxjump[ll]);
		else
			for (ll = 0; ll<NV_PH15; ll++)
				thetaH1[ll] -= min(max(NR_Numerator[ll] / NR_Denominator[ll], -maxjump[ll]), maxjump[ll]);

		// Any other constraints
		thetaH1[3] = max(thetaH1[3], 1.0f);
		thetaH1[4] = max(thetaH1[4], 0.1f);

		thetaH1[5] = max(thetaH1[5], 0.5f);
		thetaH1[5] = min(thetaH1[5], szXY / 2.0f);

		thetaH1[6] = max(thetaH1[6], 0.5f);
		thetaH1[6] = min(thetaH1[6], szZ / 2.0f);
	}

	//Estimate background model   
	thetaH0[0] = 0.0;
	for (hh = 0; hh < szZ; hh++) for (ii = 0; ii<szXY; ii++) for (jj = 0; jj<szXY; jj++) {
		//TODO add oo estimation!!
		thetaH0[0] += s_data[hh*szXY*szXY + szXY*ii + jj];
	}
	thetaH0[0] = thetaH0[0] / pow((float)szXY, 2) / szZ;
	
	// Calculating the CRLB and LogLikelihood
	Div = 0.0;
	for (hh = 0; hh < szZ; hh++) for (ii = 0; ii<szXY; ii++) for (jj = 0; jj<szXY; jj++) {
		PSFx = kernel_IntGauss1D(ii, thetaH1[0], thetaH1[5]);
		PSFy = kernel_IntGauss1D(jj, thetaH1[1], thetaH1[5]);
		PSFz = kernel_IntGauss1D(hh, thetaH1[2], thetaH1[6]);

		model = thetaH1[4] + max(thetaH1[3], thetaH1[4])*PSFx*PSFy*PSFz;
		data = s_data[hh*szXY*szXY + szXY*ii + jj];

		//calculating derivatives
		kernel_DerivativeIntGaussPSF1D(ii, thetaH1[0], thetaH1[5], thetaH1[3], PSFy, PSFz, &dudt[0], NULL);
		kernel_DerivativeIntGaussPSF1D(jj, thetaH1[1], thetaH1[5], thetaH1[3], PSFx, PSFz, &dudt[1], NULL);
		kernel_DerivativeIntGaussPSF1D(hh, thetaH1[2], thetaH1[6], thetaH1[3], PSFx, PSFy, &dudt[2], NULL);
		dudt[3] = PSFx*PSFy*PSFz;
		dudt[4] = 1.0f;

		kernel_DerivativeIntGaussPSF3DSigma(ii, jj, kk, thetaH1[0], thetaH1[1], thetaH1[2],
			thetaH1[5], thetaH1[6], thetaH1[3], PSFx, PSFy, PSFz, &dudt[5], &dudt[6], NULL, NULL);


		//Building the Fisher Information Matrix
		for (kk = 0; kk<NV_PH15; kk++)for (ll = kk; ll<NV_PH15; ll++){
			M[kk*NV_PH15 + ll] += dudt[ll] * dudt[kk] / model;
			M[ll*NV_PH15 + kk] = M[kk*NV_PH15 + ll];
		}

		//LogLikelyhood
		logModel = model / (thetaH0[0] + 1e-5);
		if (logModel>0 && data > 0)
			Div += 2 * (data*log(logModel + 1e-5) - model + thetaH0[0]);
	}

	//// Matrix inverse (CRLB=F^-1) and output assigments
	kernel_MatInvN(M, Minv, Diag, NV_PH15);
	kernel_CalcLLRProp(Diag[3], thetaH1[3], Div, LR);

	for (kk = 0; kk<NV_PH15; kk++)
		d_Parameters[kk + (NV_PH15 + NV_PH05)*(BlockSize*bx + tx)] = thetaH1[kk];
	for (kk = 0; kk<NV_PH05; kk++)
		d_Parameters[(NV_PH15 + kk) + (NV_PH15 + NV_PH05)*(BlockSize*bx + tx)] = thetaH0[kk];
	for (kk = 0; kk<Q_P5; kk++)
		d_LogLikelihood[kk + Q_P5 * (BlockSize*bx + tx)] = LR[kk];
	for (kk = 0; kk < NV_PH15; kk++)
		d_CRLBs[kk + NV_PH15*(BlockSize*bx + tx)] = Diag[kk];

	return;
}


//*******************************************************************************************
__global__ void kernel_MLEFitSigmaXX(const float *d_data, float PSFSigma, int sz, int iterations,
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood, int Nfits){
	/*!
	* \brief basic MLE fitting kernel.  No additional parameters are computed.
	* \param d_data array of subregions to fit copied to GPU
	* \param PSFSigma sigma of the point spread function
	* \param sz nxn size of the subregion to fit
	* \param iterations number of iterations for solution to converge
	* \param d_Parameters array of fitting parameters to return for each subregion
	* \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	* \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	* \param Nfits number of subregions to fit
	*/

	//__shared__ float s_data[MEM];
	float M[NV_PS2H1*NV_PS2H1], Diag[NV_PS2H1], Minv[NV_PS2H1*NV_PS2H1];
	int tx = threadIdx.x;
	int bx = blockIdx.x;
	int BlockSize = blockDim.x;
	int ii, jj, kk, ll;
	float model, cf, df, data;
	float Div;
	float PSFy, PSFx;
	float dudt[NV_PS2H1];
	float d2udt2[NV_PS2H1];
	float NR_Numerator[NV_PS2H1], NR_Denominator[NV_PS2H1];
	float thetaH1[NV_PS2H1];
	float thetaH0[NV_PS2H0];
	float maxjump[NV_PS2H1] = { 1e0f, 1e0f, 1e2f, 2e0f, 1e-1f, 1e-1f };
	float g[NV_PS2H1] = { 1.0f, 1.0f, 0.5f, 1.0f, 1.0f, 1.0f };
	float Nmax;
	float logModel;
	float LR[Q_PS2];

	//Prevent read/write past end of array
	if ((bx*BlockSize + tx) >= Nfits) return;

	memset(M, 0, NV_PSH1*NV_PSH1*sizeof(float));
	memset(Minv, 0, NV_PSH1*NV_PSH1*sizeof(float));

	//load data
	const float *s_data = d_data + (sz*sz*bx*BlockSize + sz*sz*tx);

	//initial values
	kernel_CenterofMass2D(sz, s_data, &thetaH1[0], &thetaH1[1], 0, 0);
	kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &thetaH1[3]);
	thetaH1[2] = max(0.0f, (Nmax - thetaH1[3]) * 2 * pi*PSFSigma*PSFSigma);
	thetaH1[4] = PSFSigma;
	thetaH1[5] = PSFSigma;
	for (kk = 0; kk<iterations; kk++) {//main iterative loop

		//initialize

		memset(NR_Numerator, 0, NV_PS2H1*sizeof(float));
		memset(NR_Denominator, 0, NV_PS2H1*sizeof(float));

		for (ii = 0; ii<sz; ii++) for (jj = 0; jj<sz; jj++) {
			PSFx = kernel_IntGauss1D(ii, thetaH1[0], thetaH1[4]);
			PSFy = kernel_IntGauss1D(jj, thetaH1[1], thetaH1[5]);

			model = thetaH1[3] + thetaH1[2] * PSFx*PSFy;
			data = s_data[sz*jj + ii];

			//calculating derivatives
			kernel_DerivativeIntGauss1D(ii, thetaH1[0], thetaH1[4], thetaH1[2], PSFy, &dudt[0], &d2udt2[0]);
			kernel_DerivativeIntGauss1D(jj, thetaH1[1], thetaH1[5], thetaH1[2], PSFx, &dudt[1], &d2udt2[1]);
			kernel_DerivativeIntGauss1DSigma(ii, thetaH1[0], thetaH1[4], thetaH1[2], PSFy, &dudt[4], &d2udt2[4]);
			kernel_DerivativeIntGauss1DSigma(jj, thetaH1[1], thetaH1[5], thetaH1[2], PSFx, &dudt[5], &d2udt2[5]);
			dudt[2] = PSFx*PSFy;
			d2udt2[2] = 0.0f;
			dudt[3] = 1.0f;
			d2udt2[3] = 0.0f;

			cf = 0.0f;
			df = 0.0f;
			if (model>10e-3f) cf = data / model - 1;
			if (model>10e-3f) df = data / pow(model, 2);
			cf = min(cf, 10e4f);
			df = min(df, 10e4f);

			for (ll = 0; ll<NV_PS2H1; ll++){
				NR_Numerator[ll] += dudt[ll] * cf;
				NR_Denominator[ll] += d2udt2[ll] * cf - pow(dudt[ll], 2)*df;
			}
		}

		// The update
		for (ll = 0; ll<NV_PS2H1; ll++)
			thetaH1[ll] -= g[ll] * min(max(NR_Numerator[ll] / NR_Denominator[ll], -maxjump[ll]), maxjump[ll]);

		// Any other constraints
		thetaH1[2] = max(thetaH1[2], 1.0f);
		thetaH1[3] = max(thetaH1[3], 0.01f);
		thetaH1[4] = max(thetaH1[4], PSFSigma / 10.0f);
		thetaH1[5] = max(thetaH1[5], PSFSigma / 10.0f);
	}

	thetaH0[0] = 0.0;
	for (ii = 0; ii<sz; ii++) for (jj = 0; jj<sz; jj++) {
		thetaH0[0] += s_data[sz*jj + ii];
	}
	thetaH0[0] = thetaH0[0] / pow((float)sz, 2);


	// Calculating the CRLB and LogLikelihood
	Div = 0.0f;
	for (ii = 0; ii<sz; ii++) for (jj = 0; jj<sz; jj++) {

		PSFx = kernel_IntGauss1D(ii, thetaH1[0], thetaH1[4]);
		PSFy = kernel_IntGauss1D(jj, thetaH1[1], thetaH1[5]);

		model = thetaH1[3] + thetaH1[2] * PSFx*PSFy;
		data = s_data[sz*jj + ii];

		//calculating derivatives
		kernel_DerivativeIntGauss1D(ii, thetaH1[0], thetaH1[4], thetaH1[2], PSFy, &dudt[0], NULL);
		kernel_DerivativeIntGauss1D(jj, thetaH1[1], thetaH1[5], thetaH1[2], PSFx, &dudt[1], NULL);
		kernel_DerivativeIntGauss1DSigma(ii, thetaH1[0], thetaH1[4], thetaH1[2], PSFy, &dudt[4], NULL);
		kernel_DerivativeIntGauss1DSigma(jj, thetaH1[1], thetaH1[5], thetaH1[2], PSFx, &dudt[5], NULL);
		dudt[2] = PSFx*PSFy;
		dudt[3] = 1.0f;

		//Building the Fisher Information Matrix
		for (kk = 0; kk<NV_PS2H1; kk++)for (ll = kk; ll<NV_PS2H1; ll++){
			M[kk*NV_PS2H1 + ll] += dudt[ll] * dudt[kk] / model;
			M[ll*NV_PS2H1 + kk] = M[kk*NV_PS2H1 + ll];
		}

		//LogLikelyhood
		logModel = model / (thetaH0[0] + 1e-5);
		if (logModel>0 && data > 0)
			Div += 2 * (data*log(logModel + 1e-5) - model + thetaH0[0]);
	}

	// Matrix inverse (CRLB=F^-1) and output assigments
	kernel_MatInvN(M, Minv, Diag, NV_PS2H1);
	kernel_CalcLLRProp(Diag[2], thetaH1[2], Div, LR);

	//write to global arrays
	//for (kk = 0; kk<NV_PS2H1; kk++) 
	//	d_Parameters[Nfits*kk + (NV_PSH1 + NV_PSH0)*BlockSize*bx + tx] = thetaH1[kk];
	//for (kk = 0; kk<NV_PS2H0; kk++) 
	//	d_Parameters[Nfits*(NV_PS2H1 + kk) + BlockSize*bx + tx] = thetaH0[kk];
	//for (kk = 0; kk<Q_PS2; kk++) 
	//	d_LogLikelihood[Nfits*kk + BlockSize*bx + tx] = LR[kk];
	//for (kk = 0; kk<NV_PS2H1; kk++) 
	//	d_CRLBs[Nfits*kk + BlockSize*bx + tx] = Diag[kk];

	for (kk = 0; kk<NV_PS2H1; kk++)
		d_Parameters[kk + (NV_PS2H1 + NV_PS2H0)*(BlockSize*bx + tx)] = thetaH1[kk];
	for (kk = 0; kk<NV_PS2H0; kk++)
		d_Parameters[(NV_PS2H1 + kk) + (NV_PS2H1 + NV_PS2H0)*(BlockSize*bx + tx)] = thetaH0[kk];
	for (kk = 0; kk<Q_PS2; kk++)
		d_LogLikelihood[kk + Q_PS2 * (BlockSize*bx + tx)] = LR[kk];
	for (kk = 0; kk<NV_PS2H1; kk++)
		d_CRLBs[kk + NV_PS2H1*(BlockSize*bx + tx)] = Diag[kk];
	return;
}
