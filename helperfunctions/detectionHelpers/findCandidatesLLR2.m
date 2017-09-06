function [ multiFrameRatioSampling, PFA,imageCrop,adj_p,PD ] = findCandidatesLLR2(image_out, PSFSigma, GLLRSignificance, PFA )
%createDStormSim Creates DStorm Simlation
% SYNOPSIS:
%   [ multiFrame, multiFrameRatioSampling ] = findCandidatesLLR(image_out, PSFSigma, ratioSampling, GLLRSignificance )
% 
% PARAMETERS:
%   image_out: Noised data
%   PSFSigma: Sigma of the PSF
%   ratioSampling: Number of fits is proportional to ratio sampling (1/ratioSampling)^2 is the number of fits per pixel.
%   GLLRSignificance: False discovery rate (FDR) controlled significance.
% 
% OUTPUTS:
%   multiFrame: resampled image average decision per pixel
%   multiFrameRatioSampling: decision per pixel
% 
%   SEE ALSO:
%       dipDStormGetDefaultParams, createDStormSim


a = 3*(2*PSFSigma+1)/2;
whole = floor(a);
part = a-whole;
ratioSampling=1;
if part == 0;
    Npixels=3*(2*PSFSigma+1)+1;
else
    Npixels=3*(2*PSFSigma+1);    
end

szx = ceil(Npixels/2):ratioSampling:size(image_out,2)-floor(Npixels/2);
szy = ceil(Npixels/2):ratioSampling:size(image_out,1)-floor(Npixels/2);

if nargin < 5 || isempty(PFA)
    [Ym,Xm,Zm] = meshgrid(szy-1,szx-1,0:size(image_out,3)-1);
    XX = reshape(Xm,[ numel(Xm) 1]);
    YY = reshape(Ym,[ numel(Ym) 1]);
    ZZ = reshape(Zm,[ numel(Zm) 1]);
       
    %initialize variables to retain raw fitting results
    info.Raw0=[];
    info.CRLB_STD=[];
    info.LL=[];
    info.info = []; 

    st =1;
    en = length(XX);
    iterations=8;
    disp('cMakeSubregions')
    [ stack ] = cMakeSubregions(XX(st:en),YY(st:en),ZZ(st:en),Npixels,single(image_out));
    disp('gpuGaussMLEv3')
    [P, C, L]=gpuGaussMLEv3(permute(single(stack),[2 1 3]),single(PSFSigma),iterations,0,false); 
    info.LL=[info.LL; L']
    info.Raw0=[info.Raw0; P'];
    info.CRLB_STD=[info.CRLB_STD; C'];    
    PFA=info.LL(:,3);
    PD = L(4,:);
end

imageCrop = double(image_out(szx,szy,:));
[ h, ~, adj_p ]=fdr_bh(PFA,GLLRSignificance,'dep','no');
adj_p = reshape(adj_p',[size(szx,2) size(szy,2) size(image_out,3)]);
multiFrameRatioSampling = logical(reshape(h',[size(szx,2) size(szy,2) size(image_out,3)]));
end

