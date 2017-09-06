% Code for the publication:
%  C.S. Smith, S. Stallinga, K.A. Lidke, B. Rieger, D. Grünwald, 
%  Probability-based particle detection that enables threshold-free and robust in vivo single molecule tracking,
%  Molecular Biology of the Cell, 2015. mbc.E15-06-0448
%
% Copyright 2015, University of Massachusetts, United States of America
%
% Carlas Smith, October 2015
%
% This code requires the toolbox DIPimage and CUDA 6.5 to be installed
% The windows watchdog that will stop the GPU execution if not disabled! 
% (TdrLevel=0,TdrDelay=180; https://msdn.microsoft.com/en-us/library/windows/hardware/ff569918(v=vs.85).aspx) 
% (http://www.diplib.org/download; https://developer.nvidia.com/cuda-toolkit-65)

%Bugfix:  5/18/2016, support of non-square data
%Feature: 6/7/2017, change of colour coding and probablity scale
%Bugfix: 9/4/2017, reduced memory uses

% close all
clearvars
addpath(genpath('helperfunctions'))

%% Load data

fileStructure(1).Geller = './data/grid_082015-1.tif';
fileStructure(1).Dark = './data/dark-1.tif';
fileStructure(1).Data = './data/merged3-1.tif';

%% camera calibration
%out(2) the 1/slope of the fit in [photons per ADU]: conversion ADU =
%detected #photons * slope !

a = bfopen( [ fileStructure(1).Geller]);
geller = double(cell2mat(permute(a{1}(1:min(100,end),1),[3 2 1])));

a = bfopen( [ fileStructure(1).Dark] );
bg = double(cell2mat(permute(a{1}(1:min(100,end),1),[3 2 1])));

%%
PSFSigma=1.39; %size of PSF
nFrames = 20; %maximum number of frames to be processed

a = bfopen( [ fileStructure(1).Data] );
data = double(cell2mat(permute(a{1}(1:min(end,nFrames),1),[3 2 1])));

%% Camera callibration option 1
% Tradional method need aditional data in which offset and gain are
% constant
out = cal_readnoise(geller, bg);
dataTemp=(data-repmat(mean(bg,3),[1 1 size(data,3)]))*out(2);

% %% Camera callibration option 2
% 
% for i=1:10
%     [g(i),o(i)] = pcfo(data(1:450,1:512,i));
% end
% 
% dataTemp = (data-mean(o))./mean(g);

%% GLRT Based Detection
dataCam1 = dataTemp(1:450,1:512,1:min(end,nFrames));
clear dataTemp


[coordsCam1,detParCam1,cutProcessCam1] = LLRMapv2(dataCam1,PSFSigma,[],-2);

% NOTE:
% % If you do not want to calculate cluster features you can specify: minPixels =[]
% [coordsCam1,detParCam1,cutProcessCam1] = LLRMapv2(dataCam1,PSFSigma,[]);
% 
% % If you want to compute the GLRT for all pixels without rejecting large areas with background compReduction=[]
% [coordsCam1,detParCam1,cutProcessCam1] = LLRMapv2(dataCam1,PSFSigma,0,[]);


%% Pre-Filter Detection Clusters Cam1
paramsPreFilterFits = getDefaultParamsPreFilterFits;
% paramsPreFilterFits = 
%     circularityMax: 3
%     circularityMin: 0.5000
%             pH1Max: 1
%             pH1Min: 0
%       minPixelDist: 0
%     clusterSizeMin: 0
%     clusterSizeMax: 50

% Set probablity of false detection to 0.05%

paramsPreFilterFits.pH1Min = 0.95;
[ maskPreFiltCam1 ] =  preFilterFits(coordsCam1,detParCam1,paramsPreFilterFits);

% Plot pre-filtered spots
optionsLLR = dipSubLoc2DSetOptions;
optionsLLR.BoxCenters = [coordsCam1(maskPreFiltCam1,2) coordsCam1(maskPreFiltCam1,1) coordsCam1(maskPreFiltCam1,3)];
optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
C=jet(256); 
% index = double(floor(255.*(detParCam1.pH1(1,maskPreFiltCam1)))+1);
index = double(floor(255.*(1-(1-detParCam1.pH1(maskPreFiltCam1))./(1-paramsPreFilterFits.pH1Min)))+1);
optionsLLR.BoxColor = C(index,:);
optionsLLR.plotBoxes = 1;
optionsLLR.im = dip_image(cutProcessCam1)';

dipSubLoc2D(optionsLLR);

%% MLE Fit Intensities Cam1
paramsFit = getDefaultParamsFit;
% paramsFit       
%     FitSigma: 0
%     Iterations: 10
%     MaxCudaFits: 100000
%     PSFSigma: 1.3900
%     BoxSize: 11.3400
coodsUnCut=round(coordsCam1+(1.5*(2*PSFSigma+1)-0.5).*[ones(size(coordsCam1,1),1) ones(size(coordsCam1,1),1) zeros(size(coordsCam1,1),1)]);
paramsFit.FitSigma=true;

% Fit detections
[ rawFitResultsCam1 ] = fitBoxCenters( single(squeeze(dataCam1)),[coodsUnCut(:,2) coodsUnCut(:,1) coodsUnCut(:,3)],paramsFit);

% Display 3D matrix for which each 2D slice is fitted with Apixelated-gaussian
dipshow(rawFitResultsCam1.ROIStack,'lin')

%%
paramsFilterFits = getDefaultParamsFilterFits;
% paramsFilterFits
%       MinCRLBSTD: 0
%       MinPhotons: 0
%            MinBg: 0
%       MaxCRLBSTD: Inf
%        MaxPFAValue: 0.05
%       MaxPhotons: Inf
%            MaxBg: Inf
%     MinPixelDist: 0

% Single Pixel Likelihood Ratio estimates can also be
% used as filtering
paramsFilterFits.MaxPFAValue=0.05;
[ maskFilt1 ] =  filterFits(rawFitResultsCam1,paramsFilterFits);

% Plot filtered spots
optionsLLR = dipSubLoc2DSetOptions;
optionsLLR.BoxCenters = [rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,1) rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,2) rawFitResultsCam1.Frame(maskPreFiltCam1&maskFilt1)];
optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
C=jet(256); 
% index = double(floor(255.*(min(rawFitResultsCam1.LL(3,:),1)))+1);
index = double(floor(255.*(min(rawFitResultsCam1.PFA(maskPreFiltCam1&maskFilt1)./paramsFilterFits.MaxPFAValue,1)))+1);
optionsLLR.BoxColor = C(index,:);
optionsLLR.plotBoxes = 1;
optionsLLR.im = dip_image(dataCam1)';
dipSubLoc2D(optionsLLR);


%% display fit results 

figure
C=jet(256); 
index = double(floor(255.*(1-(1-detParCam1.pH1(maskPreFiltCam1&maskFilt1))./(1-paramsPreFilterFits.pH1Min)))+1);
BoxColor = C(index,:);
subplot(3,2,3:4)
scatter3(rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,1),rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,2),rawFitResultsCam1.Frame(maskPreFiltCam1&maskFilt1),60,BoxColor,'marker','.')
colormap(jet);
caxis([paramsPreFilterFits.pH1Min 1]);
c = colorbar;
ylabel(c,'Detection Probability')
xlabel('x-position [pixel]')
ylabel('y-position [pixel]')
zlabel('t-frame [pixel]')
ntitle('Average Detection Probablity')

subplot(3,2,5:6)
index = double(floor(255.*(min(rawFitResultsCam1.PFA(maskPreFiltCam1&maskFilt1)./paramsFilterFits.MaxPFAValue,1)))+1);
BoxColor = C(index,:);
scatter3(rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,1),rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,2),rawFitResultsCam1.Frame(maskPreFiltCam1&maskFilt1),60,BoxColor,'marker','.')
colormap(jet);
caxis([0 paramsFilterFits.MaxPFAValue]);
c = colorbar;
ylabel(c,'False Alarm Probability')
xlabel('x-position [pixel]')
ylabel('y-position [pixel]')
zlabel('t-frame [pixel]')
ntitle('Point Estimate of False Alarm Probablity')

subplot(3,2,1)
histogram(rawFitResultsCam1.Photons,10)
ylabel('Count [#]')
xlabel('Signal [#Photons]')
subplot(3,2,2)
histogram(rawFitResultsCam1.Bg,10)
ylabel('Count [#]')
xlabel('Background [#Photons/pixel]')