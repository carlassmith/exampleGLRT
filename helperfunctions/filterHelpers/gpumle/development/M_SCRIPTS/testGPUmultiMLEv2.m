%% load simulated dataset
load('NP138g_1nMdye4-2011-2-15-17-19.mat','sequence');
imseries=sequence(:,:,0:499);
%1000 frame used in fitting. sequence have 10,000 frame in total.
%Bigger images (ie 512 by 512) would result in a higher memory comsuption
%would then potentially crush matlab (or output 0 matrix for SuperRes)
%because of overflow in host memory or GPU global memory.

%Try split frames into chucks and looping through them.
%ie: for 10000 frames, create a loop to fit 1000 each iteration for 10
%iterations and cat(x y intensity) or sum (super resolution image) all result 
%together after the loops.
%% fit with single fitting Nmax=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATTENTION: Multiemitter code isn't designed to perform the easiest way
% at single emitter fitting level although it achieves the same final result
% but you will need a lot twinking to achieve that. You should seek the ordinary way to perform single
% emitter fitting if it is desired. One of the option is a MATLAB toolbox
% package designed by Lidkelab AT UNM according to Smith and Lidke et al, nature
% method 2010. Please email Jason Byars (jbyars@salud.unm.edu) from Lidkelab for the distribution.
% You can also contact me if you find difficulty contacting members from
% Lidkelab.

data=imseries(:,:,0:499); %data is the input which should be poisson distributed for every pixel.
               %This can be achieved by doing a poisson calibration such
               %that after calibration, for each pixel, under constant
               %average intensity, mean=variance.
               
CCDGain = 25;
CCDOffset = 90;
               
data=(data-CCDOffset)/CCDGain;  % using camera info on Phys Cam

PSFSigma=1.3; %estimated PSF sigma in pixels

frames=1000;  %number of frames for each GPU run. Change to smaller number if
              %seeing empty super resolution image output.

Nave=200;     %Initial guess for intensity. Intensity is also fitted in 
              %this 1033 version unlike described in original publication 2011 May.
              %This is still an essential guess that would affect model
              %selection
              
Nmax=5;       %Maximum emitter number in fitting model. ie: 1 for single emitter fitting
              %5 for multi emitter fitting.
              
pvalue_threshold=0.01; %p_value threshold used for test goodness of fitting 
                       %using significance test
                       
resolution=100;        %in nano meters. Resolution target. Estimates with higher uncertainty value
                       %is discarded.
                       
pixelsize=106;         %in nano meters. Pixel size of the camera=real_size/magification.
                       %ie: 100 nm or 80 nm
                       
zm=10;                 %zoom factor. Avoid zoom as a variable name. 

boxsz=7;               %fitting subregion box size in pixels.

Dim1Size = size(data,1);
Dim2Size = size(data,2);

%% Apply Filters to get box centers and sub regions
sz = round(boxsz/2);
in = data;
 
thresh=(1/4)*Nave*(erf(sz/2/PSFSigma)/sz^2);
im_unif=unif(in,[sz sz 0],'rectangular')-unif(in,[2*sz 2*sz 0],'rectangular');
im_max=(im_unif>=.999*maxf(im_unif,[boxsz boxsz 0],'rectangular'))&(im_unif>thresh);

% use filtered regions to get box centers
BoxCenters=findcoord(im_max);

Dim1=BoxCenters(:,1);
Dim2=BoxCenters(:,2);
Dim3=BoxCenters(:,3);

% put data in single float format
data=single(permute(data,[2 1 3]));

% make sub regions with box center information
[ROIStack Dim1Start Dim2Start]=cMakeSubregions(Dim1,Dim2,Dim3,boxsz,data);

MaxCudaFits = 100000;

tic;
NROI=size(ROIStack,3);
%if NROI>MaxCudaFits
    Nloops=ceil(NROI/MaxCudaFits);
    Params=[];
    CRLB_STD=[];
    LL=[];
    
%    for nn=1:Nloops
%        st=(nn-1)*MaxCudaFits+1;
%        en=min(nn*MaxCudaFits,NROI);
        % call multi fit here
%        [x y t uncerx uncery uncerbg cov NLLR fitN b n]=GPUmultiMLEv2(ROIStack,PSFSigma,Nave,Nmax,pvalue_threshold,boxsz);
%        Params=cat(1,Params,[P Dim3(st:en)]);
%    end
%else
    % call multi fit here
    [Params CRLB cov NLLR fitData]=GPUmultiMLEv2(ROIStack,PSFSigma,Nave,Nmax,pvalue_threshold);
    
%end
tfit=toc;

% associate the frame time with the output region 
tframes = Dim3(fitData(:,1)+1);
% get actual positions
xreal = 0*Params(:,1);
yreal = 0*Params(:,2);
for ii = 1:size(Params,1)
    xreal(ii,1) = Params(ii,1) + Dim1Start(fitData(ii,1)+1);
    yreal(ii,1) = Params(ii,2) + Dim2Start(fitData(ii,1)+1);
end
% zoom actual positions
xzoom = single(xreal*zm);
yzoom = single(yreal*zm);

% reconstruct an image
xsize=(Dim1Size*zm);
ysize=(Dim2Size*zm);
SR_HistogramImage = cHistRecon(xsize,ysize,xzoom,yzoom,0);
SR_HistogramImage = permute(dip_image(SR_HistogramImage,'uint8'),[2 1]);
h=dipshow(SR_HistogramImage);





