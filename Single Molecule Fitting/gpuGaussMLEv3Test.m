% Run script by pressing run buton in menu! Paths will autometically be
% added

close all
clearvars
addpath(genpath('../helperfunctions'))

%% This script demonstrates the use of gpuGaussMLEv3 fittype 1

Nfits=100000;   %number of images to fit
bg=1;           %background fluorescence in photons/pixel/frame
Nphotons=100;   %expected photons/frame
Npixels=7;      %linear size of fit region in pixels. 
PSFsigma=1;     %PSF sigma in pixels
fittype=1;

%   Generate a stack of images
coords=Npixels/2-1+rand([Nfits 2])+[0*ones(Nfits,1) zeros(Nfits,1)];
[out] = finiteGaussPSFerf(Npixels,PSFsigma,Nphotons,bg,coords);

%   Corrupt with Poisson noise 
data=single(noise(out,'poisson',1)); %requires DipImage
%data = poissrnd(out); %requires statistics toolbox
%   Can look at data (DipImage)
%dipshow(permute(data,[2 1 3])); %dipimage permuted 1st two dimensions

iterations=12;
displayDebugInfo=false;
maxCudaFits = 1e6;

%   Fit and calculate speed
[P,CRLB,LL] = gpuGaussMLEv3(single(data),single(PSFsigma),iterations,fittype,displayDebugInfo,maxCudaFits); 

CRLBx=CRLB(1,:);
X=P(1,:)';


%   Report some details
s_x_found=std(X-coords(:,1));
meanstd=mean(sqrt(CRLBx));

fprintf('The standard deviation of x-position error is %g \n',s_x_found)
fprintf('The mean returned CRLB based x-position uncertainty is %g \n',meanstd)

