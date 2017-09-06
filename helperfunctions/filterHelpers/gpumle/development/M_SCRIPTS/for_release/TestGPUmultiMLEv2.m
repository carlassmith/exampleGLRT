%% This script demonstrates the use of GPUgaussMLE
clear;
close all;

Nfits=10000  %number of images to fit
bg=1;           %background fluorescence in photons/pixel/frame
Nphotons=100;   %expected photons/emitter
Npixels=9;      %linear size of fit region in pixels. 
PSFsigma=1;     %PSF sigma in pixels
Nmodel = 5;     %max number of emitters in a frame
amodel = 3;     %avg number of emitters in a frame
const = 0;      % if set to 1, every frame has avg emitters, else there is a variance to max
pvalue_threshold = 0.01;  % pvalue threshold for rejection of fits

%generate a stack of images
if const
    tNfits = Nfits*amodel;
    coords=(Npixels-1)*(rand([tNfits 2]))+[0*ones(tNfits,1) zeros(tNfits,1)];
    [out] = mfiniteGaussPSFerf(Npixels,PSFsigma,Nphotons,bg,coords);
    adjout = out(:,:,1:Nfits);
    for ii = 1:Nfits
        adjout(:,:,ii) = 0;
        adjout(:,:,ii) = sum(out(:,:,amodel*(ii-1)+1:amodel*(ii-1)+amodel),3);
        parts(amodel*(ii-1)+1:amodel*(ii-1)+amodel,1) = ii-1;
    end
    parts(:,2) = 3;
    
else
    % this section to be completed soon
    p_count = ceil(Nmodel*rand(Nfits,1));
    t_flo = sum(p_count);
    coords=(Npixels-1)*rand([t_flo 2])+[0*ones(t_flo,1) zeros(t_flo,1)];
    [out] = mfiniteGaussPSFerf(Npixels,PSFsigma,Nphotons,bg,coords);
    adjout = out(:,:,1:Nfits);
    zz = 1;
    for ii = 1:Nfits
       adjout(:,:,ii) = 0; 
       adjout(:,:,ii) = sum(out(:,:,zz:zz+p_count(ii)-1),3);
       parts(zz:zz+p_count(ii)-1,1)=ii-1;
       parts(zz:zz+p_count(ii)-1,2) = p_count(ii);
       zz = zz+p_count(ii);
    end
end

out = dip_image(adjout);
%corrupt with Poisson noise
data=noise(out,'poisson',1);

%fit and calculate speed
[P CRLB cov NLLR fitData] = GPUmultiMLEv2(permute(single(data),[2 1 3]),PSFsigma,Nphotons,Nmodel,pvalue_threshold);

CRLBx=CRLB(:,1);
CRLBy=CRLB(:,2);
X=P(:,1);
Y=P(:,2);

%fprintf('GPUgaussMLE has performed %g fits per second.\n',Nfits/t)

total = size(coords,1);

false = 0;  % false fits counter
missed = 0; % missed fits counter
matched = 0; % matched fits counter

% build our coordinate transfer place holders
locP{Nfits} = [];
cooP{Nfits} = [];
ind{Nfits} = [];
% Figure out which actual coordinates match the fit data through LAP
tic;
for ii = 1:Nfits
    locP{ii} = find(fitData(:,1) == ii-1);
    cooP{ii} = find(parts(:,1) == ii-1);
    xl = X(locP{ii});
    yl = Y(locP{ii});
    xc = coords(cooP{ii},1);
    yc = coords(cooP{ii},2);
    dm = zeros(length(locP{ii}),length(cooP{ii}));
    % build distance matrix to make proper assignments
    for jj = 1:length(locP{ii})
        for kk = 1:length(cooP{ii})
            dm(jj,kk) = (xl(jj)-xc(kk))^2+(yl(jj)-yc(kk))^2+eps;
%                       % remove possibility of large distance connections
%                       if dm(jj,kk) > 4
%                           dm(jj,kk) = 0;
%                       end
        end
    end
    % set birth and death high, its not that important, just for assignment
    b = (Npixels/3)*diag(ones(length(cooP{ii}),1));
    d = (Npixels/3)*diag(ones(length(locP{ii}),1));
    % birth and death is high, so lr can be all ones
    lr = (dm>0)';
    
    cc = [dm d; b lr];
    cc = sparse(cc);
    [links12 links21] = lap(cc);
    ind{ii} = [links12 links21];
    
    false = false + sum(links12(1:length(locP{ii}))>length(cooP{ii}));
    missed = missed + sum(links21(1:length(cooP{ii}))>length(locP{ii}));
    matched = matched + sum(links12(1:length(locP{ii}))<=length(cooP{ii}));
end

% figure out the standard deviation in x
x_diff{Nfits} = [];
for ii = 1:Nfits
    xl = X(locP{ii});
    yl = Y(locP{ii});
    xc = coords(cooP{ii},1);
    yc = coords(cooP{ii},2);
    loc_l = length(locP{ii});
    coo_l = length(cooP{ii});
    
    if loc_l == 0 || coo_l ==0
        continue;
    end
    
    temp_x_diff = 0;
    zz = 1;
    for jj = 1:loc_l
        temp_ind = ind{ii}(jj,1);
        if temp_ind > coo_l
            continue
        else
            temp_x_diff(zz) = xl(jj)-xc(temp_ind);
            zz = zz+1;
        end
    end
    x_diff{ii} = temp_x_diff;
end
x_diff_v = cell2mat(x_diff);
s_x_found=std(x_diff_v);

meanCRLBx=mean(CRLBx);
toc

fprintf('The percentage of correct fits accepted to the total is %g percent \n',100*matched/total)
fprintf('The percentage of false fits to total fits is %g percent \n',100*false/(false + matched))
fprintf('The percentage of missed localizations is %g percent \n',100*missed/total)
fprintf('The standard deviation of x-position error is %g \n',s_x_found)
fprintf('The mean returned CRLB based x-position uncertainty is %g \n',meanCRLBx)

%F=[CRLBx CRLBy CRLBn CRLBb]
%out1=out(:,:,0);
%[mean(N) mean(BG)]
%hist(N)