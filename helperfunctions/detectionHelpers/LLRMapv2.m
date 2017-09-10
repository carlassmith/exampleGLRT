function [coords,dectectionPar,cutProcess] = LLRMapv2(process,PSFSigma,minPixels,compReduction,significance,iterations,split,maxCudaFits,maxFramesPerBlock)
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
%     compReduction = 2;
%     significanceSinglePixel = 0.05;
%     split = false;
%     maxCudaFits = 1000000;
%     maxFramesPerBlock = 50;
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
        minPixels = floor((PSFSigma*1.5)^2);
    end
    if nargin < 4
        compReduction = 0;
    end
	if nargin < 5
        significance = 0.05;
    end
    if nargin < 6
        iterations = 8;
    end
    if nargin < 7
        split = true;
    end
    if nargin < 8
        maxCudaFits = 1000000;
    end
    if nargin < 9
        maxFramesPerBlock = 50;
    end

    dectectionPar=[];
    xbegin=round(1.5*(2*PSFSigma+1));
    xend=round(size(process,1)-1.5*(2*PSFSigma+1));
    ybegin=round(1.5*(2*PSFSigma+1));
    yend=round(size(process,2)-1.5*(2*PSFSigma+1));
    szx = xbegin:1:xend;
    szy = ybegin:1:yend;
    szx=szx(1:floor(size(szx,2)/4)*4);
    szy=szy(1:floor(size(szy,2)/4)*4);
    
    Nloops=ceil(size(process,3)/maxFramesPerBlock); %number of loops
    pfa =ones([length(szy) length(szx) size(process,3)]);
    hh =false([length(szy) length(szx) size(process,3)]);
    disp('Performing H1/H0 tests...')
    for nn = 1:Nloops
        st=(nn-1)*maxFramesPerBlock+1;
        en=min(nn*maxFramesPerBlock,size(process,3));
        idx = st:en;
        [Ym,Xm,Zm] = meshgrid(szx,szy,idx);
        
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
            cutProcess = process(szx,szy);
            flags =[1 1];
        elseif ndims(process) == 3
            cutProcess = process(szx,szy,idx);
            flags =[1 1 0];
        else
            error('Only 2D and 3D image stacks are supported')
        end


        cutProcess =max(cutProcess-double(smooth(cutProcess,flags*sqrt(n_smooth)*9)),0);

        locIm = ones(size(cutProcess));
        if ~isempty(compReduction)
            kernel(1).filter=g{1};
            kernel(1).origin = 3;
            kernel(2).filter=g{1};
            kernel(2).origin = 1; %3

            if ndims(cutProcess) == 3
                kernel(3).filter=1;    
            end
            V{1} = dip_separableconvolution(cutProcess,kernel,flags);
            kernel(1).filter=g{2};
            kernel(1).origin = 5;
            kernel(2).filter=g{2};
            kernel(2).origin = 1; %5
            if ndims(cutProcess) == 3
                kernel(3).filter=1;    
            end
            V{2} =dip_separableconvolution(V{1},kernel,flags);
            W{2} = double(V{1}-V{2});
            locIm =  W{2}	>mean(dip_image(W{2}),[],[1 2])+compReduction*std(dip_image(W{2}),[],[1 2]);
        end
        if ~ismatrix(cutProcess)
            idxIm = find(logical(permute(locIm,[2 1 3])));
        else
            idxIm = find(logical(permute(locIm,[2 1])));
        end
        
        [ ROIStack ] = cMakeSubregions(YY(idxIm)-1,XX(idxIm)-1,ZZ(idxIm)-1,3*(2*PSFSigma+1),single(process));
        [~, ~, LLr] = gpuGaussMLEv3(permute(single(ROIStack),[2 1 3 4]),single(PSFSigma),iterations,0,false,maxCudaFits.*ones(1,gpuDeviceCount)); 

        pfatemp =ones(size(Ym));
        hhtemp = false(size(Ym));
        
        [ ~,~,pfa_adj ]=fdr_bh(reshape(LLr(3,:),[size(LLr(3,:),2) 1]),significance,'dep','no');
        hhtemp(idxIm)=pfa_adj <= significance;
        
        if ndims(process) == 3 && size(process,3) > 1
            pfatemp(idxIm) = pfa_adj;
            pfatemp = double(min(smooth(pfatemp,PSFSigma.*[1 1 0]),1));
        elseif ndims(process) == 3 || size(process,3) == 1
            pfatemp(idxIm) = pfa_adj;
            pfatemp = double(min(smooth(pfatemp,PSFSigma),1));
        else
            error('Data has the wrong size')
        end


      pfa(:,:,idx) = pfatemp;
      hh(:,:,idx)  = hhtemp;
      disp([num2str(nn*maxFramesPerBlock) ' of ' num2str(size(process,3)) ' frames'])
    end
    coords=[];
    circularity=[];
    clusterSize=[];
    pH1=[];
    
    disp('Analysing clusters of positive pixels...')
    if ~isempty(minPixels)
        ll=newim([size(hh,2) size(hh,1) size(hh,3)],'sint32');
        for i=1:size(hh,3)
            if split && sum(sum(hh(:,:,i))) > 0
                h=squeeze(bclosing(hh(:,:,i)));
                D = -dt(h);
                D(~h) = +Inf;
                ll(:,:,i-1) = label(~dipwatershed(double(D),1).*squeeze(hh(:,:,i)),1);
            else
                ll(:,:,i-1) = label(logical(hh(:,:,i)),1);
            end
                msrResults = measure(squeeze(ll(:,:,i-1)), permute(squeeze(process(szx,szy,i)),[2 1 3]), ({'Gravity','P2A','Size'}));
                if isfield(msrResults,'Size')
                    a=msrResults.Size> minPixels;
                    subsetA = msrResults(a);
                    coords =cat(1,coords, [subsetA.Gravity'  (i-1).*ones(size( subsetA.Gravity,2),1)]);
                    circularity = cat(1,circularity,subsetA.P2A');
                    clusterSize = cat(1,clusterSize,subsetA.Size');
                    pH1 = cat(1,pH1,(min(max(2.*subsetA.Size'./((2*(PSFSigma))^2*pi),0),1)));     
                end
            if mod(i,maxFramesPerBlock) == 0
                disp(['slice ' num2str(i) ' of ' num2str(size(hh,3)) ' frames'])
            end
        end
        dectectionPar.hh=hh;
        dectectionPar.ll=ll;        
        dectectionPar.circularity = circularity;
        dectectionPar.clusterSize=clusterSize;
    else
        im_max = permute(cutProcess>=.999*maxf(cutProcess,(2*PSFSigma+1).*flags,'rectangular'),[ 2 1 3])&dip_image(hh);
        coords=findcoord(im_max);
        if ismatrix(process)
            coords=[coords zeros(size(coords,1),1)];
        end
        pH1 =1-pfa(im_max);
    end
    disp('Done!')

    dectectionPar.pfa_adj=pfa; 
    dectectionPar.pH1 = double(pH1)';
    cutProcess=process(szx,szy,:);
end