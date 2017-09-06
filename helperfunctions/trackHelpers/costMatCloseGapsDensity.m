function costMat = costMatCloseGapsDensity(movieInfo,links,st,en,costMatParams)
% COSTMATCLOSEGAPSDENSITY     cost matrix for gap closing
% 
% costMat = costMatCloseGapsDensity(movieInfo,links,st,en,costMatParams)
% 
% INPUTS
%   movieInfo - Array of size equal to the number of frames in a
%                      movie. See connectCoordsLAP for for info.
%   links - starts and ends for frame to frame links.
%   st - start indices
%   en - end indices
%   costMatParam - structure with input parameters. Use
%                  costMatCloseGapsSetOptions to set
%                  parameters important parameters are:
% 
%                   maxSearchDistPerFrame - vector with size 1 by probDim
%                                           describing the max search
%                                           distance per frame in each spatial
%                                           dimension.
%                   maxSearchDist - vector with size 1 by probDim
%                                   describing the max search distance in 
%                                   each spatial dimension 
%                   D - vector with size 1 by probDim containing expected 
%                       diffusion in each spatial dimension
%                   kon - rate for blinking on
%                   
%                   maxWvSearchDist - max search distance in the wavelength
%                                     dimension if applicable.
%                   wvJump - standard deviation of for spectral jumping.
%                   timeWindow - time window for gap closing
%                   costBD - user defined birth and death value. if zero
%                            then kon and koff are used to estimate birh
%                            and death values. 
%                   density - particle density. particles/pixels^2/lambda
% 
% OUTPUTS
%   costMat - cost matrix
%   nonlinkMarker -  
%   errFlag - 
% 
% See Also: connectCoordsLAP, costMatCloseGapsSetOptions
% 
% 
% Created by Pat Cutler August 2011 (UNM)
% Modified by Carlas Smith June 2014 (UMASS)

%% Input

%check whether correct number of input arguments was used
if nargin ~= nargin('costMatCloseGapsDensity')
    error('costMatCloseGapsDensity:narginIncorrect','costMatCloseGapsDensity: Incorrect number of input arguments!')
end

%check whether z-coordinates were input
if isfield(movieInfo,'zCoord')
    probDim = 3;
else
    probDim = 2;
end

%check for appropriate dimensions of all inputs
if ~isequal(probDim,length(costMatParams.D),...
        length(costMatParams.maxSearchDist),...
        length(costMatParams.maxSearchDistPerFrame))
    error('costMatCloseGapsDensity:inproperLengthDmaxSearchDistmaxSearchDistPerFrame',...
        ['costMatCloseGapsDensity: one of the input options ''D'' ''maxSearchDist''' ...
        'or ''maxSearchDistPerFrame'' doesn''t have the proper length of ' num2str(probDim)]);
end

%get cost matrix parameters
% D = reshape(costMatParams.D,[probDim 1]);
% kon = costMatParams.kon;
% blur = costMatParams.blur;

T = length(movieInfo);
timeWindow = costMatParams.timeWindow;

%%  build cost matrix

% This chunk here needs to go into a mex file for a speed
% boost!
births = size(links,1);
barName = sprintf('costMatCloseGapsDensity(%i tracks)',births);
disp(sprintf('costMatCloseGapsDensity(%i tracks)',births));
multiWaitbar(barName,0);

f2ftrackLength = sum(logical(links),2);

% values for sparse matrix
rr = zeros(births*floor(sqrt(births)),1);
cc = zeros(births*floor(sqrt(births)),1);
% vv = zeros(births*floor(sqrt(births)),1);
distance =zeros(births*floor(sqrt(births)),probDim);
time =zeros(births*floor(sqrt(births)),1);
pH11 =zeros(births*floor(sqrt(births)),1);
pH12 =zeros(births*floor(sqrt(births)),1);
LA1 =zeros(births*floor(sqrt(births)),probDim);
LA2 =zeros(births*floor(sqrt(births)),probDim);
spcounter = 1; % array counter

maxSearchDist = reshape(costMatParams.maxSearchDist,[probDim 1]);
maxSearchDistPerFrame = reshape(costMatParams.maxSearchDistPerFrame,[probDim 1]);

onePercent = round(births/100);
onePercent = max([onePercent 1]);
            LAenCol = [];
            LAstCol = [];

for mm=1:births
    if ~rem(mm,onePercent)
        multiWaitbar(barName,mm/births);
    end
    %get start point values
    tt2=st(mm);
    stid=links(mm,tt2);
    for nn=1:births
        
        if f2ftrackLength(nn) < costMatParams.minTrackLen ||...
                f2ftrackLength(mm) < costMatParams.minTrackLen
            continue;
        end % cut out any tracks shorter than minTrackLen
        %get end point values
        tt1=en(nn);
        if tt1>tt2;
            continue;
        end %end must precede start
        if (tt2-tt1 > timeWindow);
            continue;
        end % cut out any time gaps bigger than timeWindow
        enid=links(nn,tt1);
        % compile spatial information
        spatialEn=zeros(probDim,1); %spatial information for enid
        spatialSt=zeros(probDim,1); %spatial information for stid
        LAen=zeros(probDim,1); %localization accuracy for enid
        LAst=zeros(probDim,1); %localization accuracy for stid      
        spatialEn(1) = movieInfo(tt1).xCoord(enid,1);
        spatialSt(1) = movieInfo(tt2).xCoord(stid,1);
        spatialEn(2) = movieInfo(tt1).yCoord(enid,1);
        spatialSt(2) = movieInfo(tt2).yCoord(stid,1);
        if probDim == 3
            spatialEn(3) = movieInfo(tt1).zCoord(enid,1);
            spatialSt(3) = movieInfo(tt2).zCoord(stid,1);
        end
        
        if(size(movieInfo(tt1).xCoord,2) > 1)
            LAen(1) = movieInfo(tt1).xCoord(enid,2);
            LAst(1) = movieInfo(tt2).xCoord(stid,2);
            LAen(2) = movieInfo(tt1).yCoord(enid,2);
            LAst(2) = movieInfo(tt2).yCoord(stid,2);
            if probDim == 3
                LAen(3) = movieInfo(tt1).zCoord(enid,2);
                LAst(3) = movieInfo(tt2).zCoord(stid,2);
            end
        end
        pH1en(1) = movieInfo(tt1).pH1(enid,1);
        pH1st(1) = movieInfo(tt2).pH1(stid,1);
        
        if nn==mm;
            LAenCol(nn,:) = LAen;
            LAstCol(nn,:) = LAst;
            
            continue;
        end %can't connect to self
        
        %calculate the squared spatial distance between ends and starts
        dist = (spatialEn-spatialSt).^2;
        
        %distance allows connection?
        maxSearchDistTmp = min(maxSearchDistPerFrame*(tt2-tt1),maxSearchDist);
        if all(dist > maxSearchDistTmp.^2) %TODO all or any?
            continue;
        end

        rr(spcounter) = nn;  % identify counter locations for the rows
        cc(spcounter) = mm;  % identify counter locations for the columns

        distance(spcounter,:) = dist;
        time(spcounter) = tt2-tt1;
        LA1(spcounter,:) = LAen;
        LA2(spcounter,:) = LAst;
        pH11(spcounter) = pH1en;
        pH12(spcounter) = pH1st;
        
        spcounter = spcounter+1;
    end
end

rr(spcounter) = births;
cc(spcounter) = births;
% vv(spcounter) = 0;  % make sure matrix is full sized
% remove non indexable values
rr(rr==0) = [];
cc(cc==0) = [];

distance(length(rr)+1:end,:) = [];
time(length(rr)+1:end) = [];
LA1(length(rr)+1:end,:) = [];
LA2(length(rr)+1:end,:) = [];
pH11(length(rr)+1:end) = [];  % make sure matrix is full sized
pH12(length(rr)+1:end) = [];  % make sure matrix is full sized

time(spcounter) = 0;  % make sure matrix is full sized
distance(spcounter,:) = 0;  % make sure matrix is full sized
LA1(spcounter,:) = 0;  % make sure matrix is full sized
LA2(spcounter,:) = 0;  % make sure matrix is full sized
pH11(spcounter) = 0;  % make sure matrix is full sized
pH12(spcounter) = 0;  % make sure matrix is full sized
%%
kon = costMatParams.kon;
koff = costMatParams.koff;
D = costMatParams.D;
density = costMatParams.density;
%birth
stTmp = min(costMatParams.timeWindow,st-1);
%death
enTmp = min(costMatParams.timeWindow,T-en-1);

[ dm ] = transP(2,2,kon,koff,D,time,density,distance,probDim,LA1,LA2,rr,cc);
[ b ] = transP(1,2,kon,koff,D,stTmp,density,distance,probDim,LAenCol,LAstCol,rr,cc);
[ d ] = transP(2,1,kon,koff,D,enTmp,density,distance,probDim,LAenCol,LAenCol,rr,cc);
[ lr ] = transP(1,1,kon,koff,D,time,density,distance,probDim,LA1,LA2,rr,cc);

dm = dm-sparse(rr,cc,log(eps+pH11))-sparse(rr,cc,log(eps+pH12));
b = b-sparse(diag(diag(sparse(rr,cc,log(eps+1-pH11)))))-sparse(diag(diag(sparse(rr,cc,log(eps+pH12)))));
d = d-sparse(diag(diag(sparse(rr,cc,log(eps+pH11)))))-sparse(diag(diag(sparse(rr,cc,log(eps+1-pH12)))));
lr = lr-sparse(rr,cc,log(eps+1-pH11))-sparse(rr,cc,log(eps+1-pH12));

lr(isinf(lr)) = 0;
b(isinf(b)) = 0;
d(isinf(d)) = 0;
dm(isinf(dm)) = 0;    
costMat=[sparse(double(dm)) d;b lr];
costMat(isinf(costMat)) = 0;



%% Running a cost matrix check to look for questionable connections
% Any connection that doesn't hold truth after a slight perturbation is
% disregarded as too uncertain to assure a good connection

[links12, ~] = lap(real(costMat),0,0,0);
mask = ones(size(dm));

tempcMat = costMat;
for ii = 1:length(links12)
    if ii <= size(mask,1)
        if links12(ii) <= size(mask,2)
            tempcMat(ii,links12(ii)) = tempcMat(ii,links12(ii)) + costMatParams.blur; %0.69;
        end
    end
end

[tlinks12, ~] = lap(real(tempcMat),0,0,0);

if nnz(tlinks12 - links12) ~= 0
    remlinks = find(tlinks12 - links12);
    
    for jj = 1:length(remlinks)
       zz = remlinks(jj);
       if zz <= size(mask,1)
           mask(zz,:) = 0;
           if tlinks12(zz) <= size(mask,2)
               mask(:, tlinks12(zz)) = 0;     
           end
           if links12(zz) <= size(mask,2)
               mask(:, links12(zz)) = 0;
           end
       end
    end
    
    %build cost matrix
    costMat=[sparse(double(dm.*mask)) d;b lr.*mask'];
    costMat(isinf(costMat)) = 0;
    
end

multiWaitbar(barName,'close');



