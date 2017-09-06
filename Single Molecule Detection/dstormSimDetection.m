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

clearvars
close all
addpath(genpath('helperfunctions'))
addpath(genpath('../helperfunctions/'))

%%
options = dipDStormGetDefaultParams;
options.camera_params.fov = [64 64]; % number of pixels in field of vieuw

% Acquisition parameters
options.acquisition_params.timeframes = 10; % Duration of the acquisition
PSFSigma = 0.3*options.sample_params.lambda/options.optics_params.NA/options.camera_params.ps; %size of gaussian PSF

%create a siemens star
n_arms=6; % number of arms of siemens star
object_im = (cos(phiphi(options.camera_params.fov)*n_arms)>zeros(options.camera_params.fov));
minPix =  1.5*(2*PSFSigma+1)-0.5;
maxPix = options.camera_params.fov(1)-1.5*(2*PSFSigma+1)-0.5;
mask = newim(options.camera_params.fov);
mask(ceil(minPix)-1:floor(maxPix),ceil(minPix)-1:floor(maxPix)) = 1.0;
options.object_im = object_im & mask;

% density of single molcules
rho = 750;

% Fluorophore properties
% Emitter switching rate constants
options.sample_params.k_off = 1;
options.sample_params.N=sum(options.object_im )*(options.camera_params.ps/1000)^2*rho;
options.sample_params.k_on = options.sample_params.k_off/(5*rho);
        
% Create simulated data
[ imgNoiseFree, emitterStatesNoiseFree, GTImage, numOfTrueP] = createDStormSim( options );    

options.sample_params.photons = 500;    % Expected collected photons per full frame from single molecule
options.sample_params.bg = 10;          % Background [photons/pixel]    

[image_out, emitter_states] = addNoiseDStormSim(imgNoiseFree,emitterStatesNoiseFree,options);

% show simulation
dipshow(imgNoiseFree,'lin')
dipshow(image_out,'lin')

% Do GLRT detection on all frames
[coordsCam1,detParCam1,cutProcessCam1] = LLRMapv2(image_out,PSFSigma,0,[]);

% NOTE:
% % If you do not want to calculate cluster features you can specify: minPixels =[]
% [coordsCam1,detParCam1,cutProcessCam1] = LLRMapv2(dataCam1,PSFSigma,[]);
% 
% % If you want to compute the GLRT for all pixels without rejecting large areas with background compReduction=[]
% [coordsCam1,detParCam1,cutProcessCam1] = LLRMapv2(dataCam1,PSFSigma,0,[]);

% show detection result
dipshow(extend(detParCam1.hh,size(image_out)),'lin')

%% Display detection on noisy data
optionsLLR = dipSubLoc2DSetOptions;
optionsLLR.BoxCenters = [coordsCam1(:,1) coordsCam1(:,2) coordsCam1(:,3)];
optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
C=jet(256); 
index = double(floor(255.*(detParCam1.pH1))+1);
optionsLLR.BoxColor = C(index,:);
optionsLLR.plotBoxes = 1;
optionsLLR.im = dip_image(cutProcessCam1,'uint8')';

h = dipSubLoc2D(optionsLLR); 

%% Display detection on ground truth 
optionsLLR = dipSubLoc2DSetOptions;
optionsLLR.BoxCenters = [coordsCam1(:,1) coordsCam1(:,2) coordsCam1(:,3)];
optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
C=jet(256); 
index = double(floor(255.*(detParCam1.pH1(:)))+1);
optionsLLR.BoxColor = C(index,:);
optionsLLR.plotBoxes = 1;
optionsLLR.im = cut(imgNoiseFree,size(cutProcessCam1));

h = dipSubLoc2D(optionsLLR);  

