% Run script by pressing run buton in menu! Paths will autometically be
% added if the the necessary example data isdownloaded
% (www.umassmed.edu/grunwaldlab/) and places in the empty folder: data

% Code for the publication:
%  C.S. Smith, S. Stallinga, K.A. Lidke, B. Rieger, D. Grünwald, 
%  Probability-based particle detection that enables threshold-free and robust in vivo single molecule tracking,
%  Molecular Biology of the Cell, 2015 accepted. mbc.E15-06-0448;
%  L.J. van Vliet, D. Sudar and I.T. Young in Cell Biology volume III, J.E. Celis (eds.)
%    Digital Fluorescence Imaging Using Cooled CCD Array Cameras, pp.109-120, 1998
%
% Copyright 2015, University of Massachusetts, USA
%
% Carlas Smith, October 2015
%
% This code requires the toolbox DIPimage to be installed
% (http://www.diplib.org/download)

close all
clearvars
addpath(genpath('../helperfunctions'))

%% This script demonstrates how to calibrate an EMCCD gain and subtract dark counts.


fileStructure(1).Geller = '../data/grid_082015-1.tif';
fileStructure(1).Dark = '../data/dark-1.tif';
fileStructure(1).Data = '../data/merged3-1.tif';

%% intensity calibration

a = bfopen( [ fileStructure(1).Geller]);
geller = double(cell2mat(permute(a{1}(1:min(100,end),1),[3 2 1])));

a = bfopen( [  fileStructure(1).Dark] );
bg = double(cell2mat(permute(a{1}(1:min(100,end),1),[3 2 1])));

% estimate EMCCD gain
out = cal_readnoise(geller, bg);

% show the data used for calibration
dipshow(bg,'lin')
dipshow(geller,'lin')

% Caluclate callibrated image. Subtract counts and devide by gain.
%out(2) the 1/slope of the fit in [photons per ADU]: conversion ADU =
%detected #photons * slope !

nFrames=20; %Number of frames to use for data processing.
a = bfopen( [ fileStructure(1).Data] );
data = double(cell2mat(permute(a{1}(1:min(end,nFrames),1),[3 2 1])));
data=(data-repmat(mean(bg,3),[1 1 size(data,3)]))*out(2);

dipshow(data,'lin')

