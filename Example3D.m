% Code for the publication:
%  C.S. Smith, S. Stallinga, K.A. Lidke, B. Rieger, D. Grünwald, 
%  Probability-based particle detection that enables threshold-free and robust in vivo single molecule tracking,
%  Molecular Biology of the Cell, 2015 accepted. mbc.E15-06-0448
%
% Copyright 2017, University of Oxford, United Kingdom
%
% Carlas Smith, June 2017
%
% This code requires the toolbox DIPimage and CUDA 6.5 to be installed
% The windows watchdog that will stop the GPU execution if not disabled! 
% (TdrLevel=0,TdrDelay=180; https://msdn.microsoft.com/en-us/library/windows/hardware/ff569918(v=vs.85).aspx) 
% (http://www.diplib.org/download; https://developer.nvidia.com/cuda-toolkit-65)

close all
clearvars
addpath(genpath('helperfunctions'))

%%
PSFSigma=1.39; %nyquist size of PSF sigma in pixels
numberOfChannels = 4;

fileName = './data/20170326_dlgYFP_YFPATTO647N_dlg570_p3s5r.ome.tiff';
a = bfopen(fileName);
img = double(cell2mat(permute(a{1}(:,1),[3 2 1])));
img = reshape(img,[size(img,1) size(img,2)   size(img,3)/numberOfChannels numberOfChannels]);

data670 = img(:,:,:,1);

%% Camera callibration option 2
% New method uses a single image to perform callibration
for i=1:10
    [g(i),o(i)] = pcfo(data670(1:450,1:512,i));
end

data670 = (data670-mean(o))./mean(g);


channel=2;  
numberOfChannels = 4;

fileName = './data/20170326_dlgYFP_YFPATTO647N_dlg570_p3s5r.ome.tiff';
a = bfopen(fileName);
img = double(cell2mat(permute(a{1}(:,1),[3 2 1])));
img = reshape(img,[size(img,1) size(img,2)   size(img,3)/numberOfChannels numberOfChannels]);

data670 = img(:,:,:,1);

fileName = './data/20170326_640_Calibration1.ome.tiff';
a = bfopen(fileName);
Calibration640 = double(cell2mat(permute(a{1}(:,1),[3 2 1])));

fileName = './data/20170326_640_darkCalibration1.ome.tiff';
a = bfopen(fileName);
darkCalibration640 = double(cell2mat(permute(a{1}(:,1),[3 2 1])));

[out640] = cal_readnoise(Calibration640,darkCalibration640);
data670 = (data670-mean(darkCalibration640,3))*out640(2);


%% GLRT Based Detection

[coordsCam1,detParCam1,cutProcessCam1] = LLRMap3D(data670,[PSFSigma 2.5*PSFSigma],[],3);

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
optionsLLR.BoxCenters = [coordsCam1(maskPreFiltCam1,2) coordsCam1(maskPreFiltCam1,1) floor(coordsCam1(maskPreFiltCam1,3))];
optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
C=jet(256); 

index = double(floor(255.*(1-(1-detParCam1.pH1(maskPreFiltCam1))./(1-paramsPreFilterFits.pH1Min)))+1);
optionsLLR.BoxColor = C(index,:);
optionsLLR.plotBoxes = 1;
optionsLLR.im = dip_image(cutProcessCam1)'; %repmat(max(cutProcessCam1,[],3),[1 1 10]))';

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
[ rawFitResultsCam1 ] = fitBoxCenters( single(squeeze(data670)),[coodsUnCut(:,2) coodsUnCut(:,1) coodsUnCut(:,3)],paramsFit);

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
paramsFilterFits.MaxPFAValue=1;
[ maskFilt1 ] =  filterFits(rawFitResultsCam1,paramsFilterFits);

% Plot filtered spots
optionsLLR = dipSubLoc2DSetOptions;
optionsLLR.BoxCenters = [rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,1) rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,2) rawFitResultsCam1.Frame(maskPreFiltCam1&maskFilt1)];
optionsLLR.BoxSize = 2*(3*PSFSigma+1).*[1 1];
C=jet(256); 

index = double(floor(255.*(min(rawFitResultsCam1.PFA(maskPreFiltCam1&maskFilt1)./paramsFilterFits.MaxPFAValue,1)))+1);
optionsLLR.BoxColor = C(index,:);
optionsLLR.plotBoxes = 1;
optionsLLR.im = dip_image(data670)';
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

% List of 3D positions x,y by MLE and z by center of mass of signification
% hypothesis tested pixels
threeDPos = [rawFitResultsCam1.Coord(maskPreFiltCam1&maskFilt1,1:2) coodsUnCut(maskPreFiltCam1&maskFilt1,3)];