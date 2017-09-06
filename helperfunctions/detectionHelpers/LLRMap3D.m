function [coords,dectectionPar,cutProcess] = LLRMap3D(process,PSFSigma,minPixels,compReduction,significance,iterations,split,maxCudaFits)
%LLRMapv2 Detects single mocules according to the GLRT framework
% SYNOPSIS:
%  [coords,dectectionPar,cutProcess] = LLRMapv2(process,PSFSigma,minPixels,
%                               compReduction,iterations,split,maxCudaFits)
% 
% PARAMETERS:
%     process: single molcule data (corrected for gain and dark counts)
% 
%     PSFSigma: sigma of diffraction limited PSF in pixels
% 
%     minPixels: minimal size of detection cluster. minPixels = [] obbits
%     calculation of detection features dectectionPar =[]
% 
%     compReduction: computational complexity reduction level in standard 
%     deviations from mean. Rejectregions with small change of single 
%     molecules, based on wavelet filtering. compReduction = [] results in 
%     computation of GLRT on all pixels.
% 
%     significance: false discovery rate 
% 
%     iterations: number of iterations for MLE algorithm.
% 
%     split: split clusters using watershed algorithm
% 
%     maxCudaFits: maximum numbers of cuda fits on single GPU device. If
%     n CUDA enables devices are ready 1xn vector is expected.
% 
% 
% DEFAULTS:
%     minPixels = (PSFSigma*1.5)^2; 
%     compReduction = 3;
%     significanceSinglePixel = 0.05;
%     split = false;
%     maxCudaFits = 1000000;
% 
% 
% OUTPUTS:
%   coords: matrix with detection coordinates.
%   dectectionPar
%         .circularity: circularity of cluster using P2A feature.
%                 .pH1: detection probability of cluster.
%         .clusterSize: cluster size
%                  .ll: labeled clusters
%                  .hh: binary image containing raw detections
%             .pfa_adj: false discovery rate corrected false positive 
%                       probailty
%             
%   cutProcess: cropped image with the actuall test pixels.

    if nargin < 3 || isempty(minPixels)
        minPixels = 0; %floor(1.5*PSFSigma(2)*(PSFSigma(1)*1.5)^2);
    end
    if nargin < 4
        compReduction = 3;
    end
	if nargin < 5
        significance = 0.05;
    end
    if nargin < 6
        iterations = 20;
    end
    if nargin < 7
        split = true;
    end
    if nargin < 8
        maxCudaFits = 1000000;
    end
    PSFSigma=PSFSigma(1);

    dectectionPar=[];
    xbegin=round(1.5*(2*PSFSigma(1)+1));
    xend=round(size(process,1)-1.5*(2*PSFSigma(1)+1));
    ybegin=round(1.5*(2*PSFSigma(1)+1));
    yend=round(size(process,2)-1.5*(2*PSFSigma(1)+1));
%     zetbegin=round(1.5*(2*PSFSigma(2)+1));
%     zetend=round(size(process,3)-1.5*(2*PSFSigma(2)+1));

    szx = xbegin:1:xend;
    szy = ybegin:1:yend;
    szzet = 1:size(process,3);
%     szzet = zetbegin:1:zetend;
    szx=szx(1:floor(size(szx,2)/4)*4);
    szy=szy(1:floor(size(szy,2)/4)*4);
%     szzet=szzet(1:floor(size(szzet,2)/4)*4);
    [Ym,Xm,Zm] = meshgrid(szx,szy,szzet);
%     scatter(reshape(Ym(:,:,1),[numel(Ym(:,:,1)) 1]),reshape(Xm(:,:,1),[numel(Xm(:,:,1)) 1]))
    XX = reshape(Xm,[ numel(Xm) 1]);
    YY = reshape(Ym,[ numel(Ym) 1]); 
    ZZ = reshape(Zm,[ numel(Zm) 1]);

    H2 = 1/16;
    H1 = 1/4;
    H0 = 3/8;
    g{1} = [H2,H1,H0,H1,H2];
    g{2} = [H2,0,H1,0,H0,0,H1,0,H2];

    n_smooth=2;
    if ismatrix(process)
%         cutProcess = double(cut(permute(dip_image(process),[2 1 3]),[size(Ym,1) size(Ym,2)]));
        cutProcess = process(szx,szy);
        flags =[1 1];
    elseif ndims(process) == 3
        cutProcess = process(szx,szy,szzet);
%         cutProcess = double(cut(permute(dip_image(process),[2 1 3]),[size(Ym,1) size(Ym,2) size(Ym,3)]));
        flags =[1 1 1];
    else
        error('Only 2D and 3D image stacks are supported')
    end
    
    
%     cutProcess =max(cutProcess-double(smooth(cutProcess,flags*sqrt(n_smooth)*9)),0);
    diskSize1=5;
    se1= strel('disk',diskSize1);
%     cutProcess = ;
% thresMult1=10;
% paramsGeneral= getDefaultParamsGeneral;
% im_max1= (out>=.999*maxf(out,2*[paramsGeneral.Psfxy paramsGeneral.Psfxy paramsGeneral.Psfz],'elliptic'));
% coordspre=findcoord(im_max1&(out>mean(dip_image(out))+std(dip_image(out))*thresMult1)); %bigger multiplication---harsh
% I570 = out(coord2image(coordspre,size(out)));
% coords570=[coordspre(:,1) coordspre(:,2) coordspre(:,3)];

    
    locIm = ones(size(cutProcess));
    if ~isempty(compReduction)
        kernel(1).filter=g{1};
        kernel(1).origin = 3;
        kernel(2).filter=g{1};
        kernel(2).origin = 1; %3
        
        if ndims(process) == 3
            kernel(3).filter=1;    
        end
        V{1} = dip_separableconvolution(imtophat(double(cutProcess),se1),kernel,flags);
        kernel(1).filter=g{2};
        kernel(1).origin = 5;
        kernel(2).filter=g{2};
        kernel(2).origin = 1; %5
        if ndims(process) == 3
            kernel(3).filter=1;    
        end
        V{2} =dip_separableconvolution(V{1},kernel,flags);
        W{2} = double(V{1}-V{2});
        locIm =  W{2}	>mean(dip_image(W{2}),[],[1 2])+compReduction*std(dip_image(W{2}),[],[1 2]);
    end
    if ~ismatrix(process)
        idxIm = find(logical(permute(locIm,[2 1 3])));
    else
        idxIm = find(logical(permute(locIm,[2 1])));
    end


%     xxs=XX(idxIm);
%     yys=YY(idxIm);    
%     imgcomp = zeros(size(process));
%     imgcomp(szx,szy,:)=locIm;
%     dipshow(max(imgcomp,[],3),'lin')
%     hold on
%     scatter(reshape(xxs(:,:,1)-1,[numel(xxs(:,:,1)) 1]),reshape(yys(:,:,1)-1,[numel(yys(:,:,1)) 1]),'marker','.')

    % paramsFit.MaxCudaFits = 2e3;
    % 
    % Psfxy=2;
    % paramsGeneral.Psfz=3;
    % paramsGeneral.Psf = single([paramsGeneral.Psfxy paramsGeneral.Psfz]);
%     Psfxy = PSFSigma(1);
%     Psfz = PSFSigma(2);
%     boxSize=single(2*[(2*Psfxy+1) (2*Psfz+1) ]);
%     sigmaPSF =  single([Psfxy Psfz]);

%     [ROIStack] = cMakeSubregionsv2(YY(idxIm)-1,XX(idxIm)-1,ZZ(idxIm)-1,ceil(boxSize(1)),single(process),ceil(boxSize(2))); %, top,left,depth
%     % XRoiStart = single(top);
%     % YRoiStart = single(left);
%     % ZRoiStart = single(depth);
% 
%     [~, ~, LLr] = gpuGaussMLEv3(permute(single(ROIStack),[2 1 3 4]),single(sigmaPSF),iterations,3,true,maxCudaFits.*ones(1,gpuDeviceCount)); 
% 
%     pfa =ones(size(Ym));
%      
%     [ ~,~,pfa_adj ]=fdr_bh(reshape(LLr(3,:),[size(LLr(3,:),2) 1]),significance,'dep','no');
%     pfa(idxIm) =double(pfa_adj);
%     pfa = min(smooth(pfa,[PSFSigma(1) PSFSigma(1) PSFSigma(2)]),1);
%     
%     coords=[];
%     circularity=[];
%     clusterSize=[];
%     pH1=[];
%     hh = false(size(Ym));
%     hh(idxIm)=pfa_adj <= significance;

    [ ROIStack ] = cMakeSubregions(YY(idxIm)-1,XX(idxIm)-1,ZZ(idxIm)-1,3*(2*PSFSigma+1),single(process));
    [~, ~, LLr] = gpuGaussMLEv3(permute(single(ROIStack),[2 1 3 4]),single(PSFSigma),iterations,0,false,maxCudaFits.*ones(1,gpuDeviceCount)); 

    pfa =ones(size(Ym));
     
    [ ~,~,pfa_adj ]=fdr_bh(reshape(LLr(3,:),[size(LLr(3,:),2) 1]),significance,'dep','no');
    if ~ismatrix(process)
        pfa(idxIm) =double(min(smooth(pfa_adj,PSFSigma),1));
    else
        pfa(idxIm) =double(min(smooth(pfa_adj,PSFSigma.*[1 1 0]),1));
    end
    coords=[];
    circularity=[];
    clusterSize=[];
    pH1=[];
    hh = false(size(Ym));
    hh(idxIm)=pfa_adj <= significance;

    if ~isempty(minPixels)
        ll=newim([size(Ym,2) size(Ym,1) size(Ym,3)],'sint32');
        if split && sum(sum(sum(hh))) > 0
            h=squeeze(bclosing(hh));
            D = -dt(h);
            D(~h) = +Inf;
            ll = label(~dipwatershed(double(D),1).*squeeze(hh),1);
        else
            ll = label(logical(hh),1);
        end
        msrResults = measure(squeeze(ll), permute(squeeze(cutProcess),[2 1 3]), ({'Gravity','P2A','Size'}));
        a=msrResults.Size> minPixels;
        subsetA = msrResults(a);
        coords =cat(1,coords, [subsetA.Gravity' ]);
        circularity = cat(1,circularity,subsetA.P2A');
        clusterSize = cat(1,clusterSize,subsetA.Size');
        pH1 = cat(1,pH1,(min(max(2.*subsetA.Size'./((2*(PSFSigma(1)))^2*pi),0),1)));     

        dectectionPar.hh=hh;
        dectectionPar.ll=ll;        
        dectectionPar.circularity = circularity;
        dectectionPar.clusterSize=clusterSize;
    else
        im_max = permute(cutProcess>=.999*maxf(cutProcess,[(2*PSFSigma(1)+1) (2*PSFSigma(1)+1) (2*PSFSigma(2)+1)],'rectangular'),[ 2 1 3])&dip_image(hh);
        coords=findcoord(im_max);
        if ismatrix(process)
            coords=[coords zeros(size(coords,1),1)];
        end
        pH1 =1-pfa(im_max); % stretch(double(),0,100,0,1);
    end

%     dipshow(hh+double(coord2image(coords,[size(hh,2) size(hh,1) size(hh,3)])),'lin')
    dectectionPar.pfa_adj=pfa; 
    dectectionPar.pH1 = double(pH1)';
end