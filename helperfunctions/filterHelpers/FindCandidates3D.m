function [cds,im_unif] = FindCandidates3D(data,PSFsigma,Photons,Pixels,ratioOfFilters)

%if input value of rationOfFilters is not given 

if nargin < 5
    ratioOfFilters=2;
end

    szx = round(Pixels(1)/2);    
    szy = round(Pixels(2)/2);    
    szz = round(Pixels(3)/2);    
    
    thresh=(1/4)*Photons*(erf(szx/2/PSFsigma(1))/szx^2)...
        *(erf(szy/2/PSFsigma(2))/szy^2)%*(erf(szz/2/PSFsigma(3))/szz^2)
    
    data(data<=0) = 1E-3;
    in=dip_image(permute(data,[2 1 3]));

    im_unif=unif(squeeze(in),[szx szy szz],'elliptic')-unif(squeeze(in),ratioOfFilters.*[szx szx szz],'elliptic');
    im_max=(im_unif>=.999*maxf(im_unif,[Pixels(1) Pixels(2) Pixels(3)],'elliptic'));
    cds=findcoord(im_max&(im_unif>thresh));
   
end