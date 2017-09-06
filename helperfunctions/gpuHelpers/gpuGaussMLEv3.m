%gpuGaussMLEv3 performs MLE fits on CUDA 6.5 enables device.
% SYNOPSIS:
%  [P,CRLB,LL] = gpuGaussMLEv3(ROIStack,PSFSigma,Iterations,fittype,display,maxCudaFits); 
%
% 
% PARAMETERS:
%     ROIStack: strack of single molecule data create
%     PSFSigma: depending on fittype scalar or vector.
%     Iterations: number of ieterations for GPU routine
%     fittype: the method of fitting is given by the fittype variable.
%       fittype=0:  Fits (Photons,Bg) under H1 and (Bg) under H0 given PSF_sigma. 
%       fittype=1:  Fits (x,y,bg,Photons) under H1 and (Bg) under H0 given PSF_sigma. 
%       fittype=2:  Fits (x,y,bg,Photons,Sigma) under H1 and (Bg) under H0 given PSF_sigma. 
%       fittype=3:  Fits (x,y,bg,Photons,z) under H1 and (Bg) under H0 given PSF_\sigma. 
%       fittype=4:  Fits (x,y,bg,Photons,Signa_x,Sigma_y) under H1 and (Bg) under H0 given PSF_\sigma. 
%       fittype=5:  Fits (x,y,z,bg,Photons,Sigma_x,Sigma_z) under H1 and (Bg) under H0 given PSF_\sigma. 
%   display: boolean to display output.
%   maxCudaFits: maximum amount of CUDA fits per GPU enables device.
% 
% 
% OUTPUTS:
%   P is the N x M matrix of found parameters 
%   CRLB is the N x M matrix of found parameters of Cramer-Rao Lower Bound calculated variances for each parameter.
%   The CRLB is now calculated internally using an LU decomposition method for inverting the Fisher information matrix.
%   The center of mass estimation is now used for all models as a starting guess for the x,y positions.  
%   This should result in fewer iterations to achieve convergence.
%   LL is the 4 x M matrix of calculated likelihood values likelihood
%   ratio, decentrality parameter, P_FA, PD, respectively.
%   where N is the number of fits and M is the number of fitted variables (i.e. 6 for fittype=4).
% 
%
% COMMENT
% This code requires the NVIDA CUDA 6.5 Driver and Toolkit to be installed. 
%  
% LITERATURE:
% Fast, single-molecule localization that achieves theoretically minimum 
% uncertainty by Carlas S. Smith, Nikolai Joseph, Bernd Rieger & Keith A.Lidke.
% 
% 
% SEE ALSO:
%   cMakeSubregions