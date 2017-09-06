% Code for the publication:
%  C.S. Smith, S. Stallinga, K.A. Lidke, B. Rieger, D. Grünwald, 
%  Probability-based particle detection that enables threshold-free and robust in vivo single molecule tracking,
%  Molecular Biology of the Cell, 2015 accepted. mbc.E15-06-0448
%
% Copyright 2015, University of Massachusetts, USA
%
% Carlas Smith, October 2015
%
% This code requires the toolbox DIPimage and CUDA 6.5 to be installed
% The windows watchdog that will stop the GPU execution if not disabled! 
% (TdrLevel=0,TdrDelay=180; https://msdn.microsoft.com/en-us/library/windows/hardware/ff569918(v=vs.85).aspx) 
% (http://www.diplib.org/download; https://developer.nvidia.com/cuda-toolkit-65)

clear all
close all
addpath(genpath('helperfunctions'))
addpath(genpath('../helperfunctions/'))

%%
bgc=20; % Background [photons/pixel]    
I0 =75; % Expected collected photons per full frame from single molecule

sigma=1.39; %PSF sigma in pixels
Npixels=ceil(3*(2*sigma+1)*2);

% Make sure we have enough pixels
part = Npixels/2-floor(Npixels/2);

if part == 0;
    Npixels=Npixels+1;
end

Nsims = 4096; %number of single molcules for simulation

% Create simulated data
coords=repmat((Npixels-1)/2,[Nsims,2]); %
M1 = finiteGaussPSFerf(Npixels,sigma,I0,bgc,coords);
dataH1=noise(M1,'poisson');
dataH0=noise(ones(Npixels,Npixels,Nsims)*bgc,'poisson');
data = cat(3,dataH1,dataH0);
significance = 0.05;

% Do GLRT single molcules detection in image
[ multiFrameRatioSampling, pfa] = findCandidatesLLR2(data, sigma, significance);

% Show results
dipshow(sum(multiFrameRatioSampling,3),'lin')
