function [ closedNucleusEdge, closedNucleus ] = getPeriphery(data,sizeNux,sizePore,nClose,nShrink,timesStdGradTreshHold,fracBackground)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    if nargin < 7
        fracBackground =10;
    end 
    
    in= mean(dip_image(data),[],[3]);
    sz = round(sizeNux/2);
    im_unif=unif(in,[sz sz 0],'rectangular')-unif(in,[2*sz 2*sz 0],'rectangular');
    im_unif = max(im_unif,max(im_unif)/fracBackground)-max(im_unif)/fracBackground;
    excOne = bdilation(im_unif > mean(im_unif) + timesStdGradTreshHold*std(im_unif),nClose);

    data(data<=0) = 1E-3;
    in= mean(dip_image(data),[],[3]);
    sz = round(sizePore/2);
    im_unif=unif(in,[sz sz 0],'rectangular')-unif(in,[2*sz 2*sz 0],'rectangular');
    excTwo = berosion(im_unif > mean(im_unif) + timesStdGradTreshHold*std(im_unif));
    labeled = squeeze((excOne & excTwo).*label(excOne));

    tempClosedNucleus= newim([size(labeled) max(labeled)],'bin');
    tempClosedNucleusEdge= newim([size(labeled) max(labeled)],'bin');
    for i=1:max(labeled)   
        tempClosedNucleus(:,:,i-1) =  berosion(bwconvhull(logical(labeled == i)),nShrink);
        tempClosedNucleusEdge(:,:,i-1)=  tempClosedNucleus(:,:,i-1)-berosion(squeeze(tempClosedNucleus(:,:,i-1)));
    end
    
    notZero = sum(tempClosedNucleus,[],[1 2]) > 0;
    closedNucleus= newim([size(labeled) sum(notZero)],'bin');
    closedNucleusEdge= newim([size(labeled) sum(notZero)],'bin');
    
    for i=1:max(labeled)   
        if notZero(i-1)
           closedNucleus(:,:,sum(notZero(0:i-1))-1) = tempClosedNucleus(:,:,i-1);
           closedNucleusEdge(:,:,sum(notZero(0:i-1))-1)=  tempClosedNucleus(:,:,i-1)-berosion(squeeze(tempClosedNucleus(:,:,i-1)));
        end
    end
	
end

