function [ tracks, out, stats, tracksStd, subRegionModel] = connectCoordsLAP(movieInfo,costMatOptions)
%CONNECTCOORDSLAP (1) links features between frames and (2) closes gaps
%
% EXAMPLE:
%   [tracksFinal,errFlag] = connectCoordsCostMat(movieInfo,costMatrices);
%
%INPUT  movieInfo    : Array of size equal to the number of frames in a
%                      movie, containing the fields:
%             .xCoord      : x-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%             .xSigma      : psf sigma for x-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation parameter estimation (zeros if not
%                            available).
%             .yCoord      : y-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%             .ySigma      : psf sigma for y-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation parameter estimation (zeros if not
%                            available).
%             .zCoord      : z-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%                            Optional. Skipped if problem is 2D. Default: zeros.
%             .zSigma      : psf sigma for z-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation parameter estimation (zeros if not
%                            available).
%             .wvCoord     : wavelength-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%                            Optional. Skipped if not defined. Default:
%                            zeros.
%             .wvSigma      : psf sigma for x-coordinates of detected features.
%                             1st column: value, 2nd column: standard
%                             deviation parameter estimation (zeros if not
%                             available).
%             .photons      : "Intensities" of detected features.
%                             1st column: values (ones if not available),
%                             2nd column: standard deviation (zeros if not
%                             available).
%             .roiIdxAll    : index for the region of interest in 
%                             from which localization originated. Relative ROIs fit
%                             for all frames. (not required)
%             .roiIdxFrame  : index for the region of interest in 
%                             from which localization originated. Relative to ROIs
%                             fit for the respective time frame (not required)
%             .model        : number of emitters fit to the model selected
%             .modelIdx     : index for specific emitter
%
%       costMatrices : struct with 2 fields indicating cost matrices and their
%                      parameters. Fields are required to be:
%                           .costMatF2F - entry supplies the cost matrix for linking
%                               between consecutive frames. Required to
%                               contain the field funcName.
%                           .costMatGC - entry supplies the cost matrix for
%                               closing gaps. Required to contain the field funcName.
%                           Note: funcName subfield is the name of function used to 
%                                 calculate cost matrix.
%
%OUTPUT 
%   out   : Structure array where each element corresponds to a
%                       track. Each element contains the following
%                       fields:
%       X: position in x dimension
%       std_x: theoretical standard deviation of position estimate x dimension
%       std_sigma_x: theoretical standard deviation for estimated
%           PSF for x dimension
%       sigma_x: estimated PSF for x dimension
%       Y: position in y dimension
%       std_y: theoretical standard deviation of position estimate y dimension
%       std_sigma_y: theoretical standard deviation for estimated
%           PSF for y dimension
%       sigma_y: estimated PSF for y dimension
%       Z: position in z dimension (dependent on zCoord input)
%       std_z: theoretical standard deviation of position estimate z dimension
%       std_sigma_z: theoretical standard deviation for estimated
%           PSF for z dimension
%       sigma_z: estimated PSF for z dimension
%       Wv: position in spectral dimension (dependent on wvCoord input)
%       std_wv: theoretical standard deviation of position estimate wv dimension
%       std_sigma_wv: theoretical standard deviation for estimated
%           PSF for wv dimension
%       sigma_wv: estimated PSF for wv dimension
%       Photons: estimated intensity in photons/frame
%       std_Photons: theoretical standard deviation for estimated
%           intensity
%       roiIdxAll: index for the region of interest in fit results
%           from which localization originated. Relative ROIs fit
%           for all frames
%       roiIdxFrame: index for the region of interest in fit results
%           from which localization originated. Relative to ROIs
%           fit for the respective time frame
%       model: number of emitters fit to the model selected
%       modelIdx: index for specific emitter
%   stats: random tracking statistics (time for key steps and TrackLinks)
%
% Dependencies: lap
%
% by Pat Cutler (UNM) August 2011

%PJC January 2012 adapted inputs/outputs for call from SPT class
% Modified by Carlas Smith May 2014 (UMASS/TU-Delft)

% Running LAP
fprintf('%s\n','Making frame-frame connections...')
t = tic;

T = length(movieInfo);
n=zeros(T,1);
births=length(movieInfo(1).xCoord);
ind=[];
cost=[];
ind{T-1} = 0;
cost{T-1} = 0;
for tt=1:T-1
    n(tt)=size(movieInfo(tt).xCoord,1);
    m=size(movieInfo(tt+1).xCoord,1);
    if n(tt) && m
         cc = costMatFrame2FrameDensity(movieInfo(tt:tt+1),costMatOptions.costMatF2F);
    else
        %2013-04-26 PJC and PKR added to deal with single observations
        if ~n(tt) && ~m
            cc = {};
        else if n(tt) && ~m
                cc = diag(ones(1,n(tt)));
            else if ~n(tt) && m
                    cc = diag(ones(1,m));
                end
            end
        end
    end
    if ~isempty(cc) && nnz(cc) ~= 0
        %submit to lap
        [links12, links21, Ass_cost] = lap(real(cc),0,0,0);
        ind{tt}=[links12 links21];
        cost{tt}=Ass_cost;
        births=births+sum(links21(1:m)>n(tt));
    end
end
n(T)=m;

links=zeros(births,T,'uint32');
%set first frame ID
links(1:n(1),1)=(1:n(1));

mt=n(1); %number of assigned tracks. This grows with each new birth.
st=zeros(births,1); %start of each track
en=zeros(births,1); %end of each track
st(1:n(1))=1;
for tt=1:T-1
    if isempty(ind{tt}) || length(ind{tt}) == 1;
        continue;
    end
    links12=ind{tt}(:,1);
    links21=ind{tt}(:,2);
    for kk=1:mt  %loop through each track
        if ~links(kk,tt);
            continue;
        end %track is dead :( (ended)
        if links12(links(kk,tt))>n(tt+1)
            links(kk,tt+1)=0; %death
            en(kk)=tt;
        else
            links(kk,tt+1)=links12(links(kk,tt)); %link to next id
        end
    end
    %add births
    for bb=1:n(tt+1)
        if links21(bb)>n(tt)    %then there was a birth (congrats!)
            %pause
            mt=mt+1;
            links(mt,tt+1)=bb;
            st(mt)=tt+1;
        end
    end
end
en(en==0)=T;

%cut out links less than minTrackLen
if isfield(costMatOptions.costMatGC,'minTrackLen')
    cm = sum(logical(links),2) < costMatOptions.costMatGC.minTrackLen;
    en(cm) = [];
    st(cm) = [];
    links(cm,:) = [];
end
births = length(en); % we only recognize births that leads to extensive tracks
fprintf('%g%s\n',births,' frame-frame tracks')
clear cm;
stats.timeFrame2Frame = toc(t);
fprintf('%s%g%s\n','Elapsed time is ',stats.timeFrame2Frame,' seconds.')

%build gap closing cost matrix
fprintf('%s\n','Making gap closing cost matrix...')
t = tic;
cc = costMatCloseGapsDensity(movieInfo,links,st,en,costMatOptions.costMatGC);
stats.timeGapCloseCostMat = toc(t);
fprintf('%s%g%s\n','Elapsed time is ',stats.timeGapCloseCostMat,' seconds.')

%% submit to lap
fprintf('%s\n','Linking gap closing cost matrix...')
t = tic;
[links12, links21, Ass_cost] = lap(real(cc),0,0,0);
clear cc
tmplinks12=links12;
tmplinks21=links21;

Ntracks=sum(links21(1:births)>births);
stats.TrackLinks=zeros(Ntracks,T,'uint32');

% Total number of tracks:
mt=1;
for nn=1:births
    id=nn;
    if ~tmplinks12(id);
        continue;
    end
    while id<=births
        stats.TrackLinks(mt,:)=stats.TrackLinks(mt,:)+links(id,:);
        id2=tmplinks12(id);
        tmplinks12(id)=0;
        id=id2;
    end
    if mt == Ntracks
        break;
    end
    mt=mt+1;
end
clear links

%check whether z-coordinates were input
if isfield(movieInfo,'zCoord')
    probDim = 3;
else
    probDim = 2;
end
%check whether wv-coordinates were input
if isfield(movieInfo,'wvCoord')
    wvDim = 1;
else
    wvDim = 0;
end

for nn=1:Ntracks
    %find all links for track nn
    TrackLinksnn = find(stats.TrackLinks(nn,:));
    count = 0; %reset counter
    %initialize output values
    out(nn).X = zeros(1,numel(TrackLinksnn));
    out(nn).std_x = zeros(1,numel(TrackLinksnn));
    out(nn).sigma_x = zeros(1,numel(TrackLinksnn));
    out(nn).std_sigma_x = zeros(1,numel(TrackLinksnn));
    out(nn).Y = zeros(1,numel(TrackLinksnn));
    out(nn).std_y = zeros(1,numel(TrackLinksnn));
    out(nn).sigma_y = zeros(1,numel(TrackLinksnn));
    out(nn).std_sigma_y = zeros(1,numel(TrackLinksnn));
    if probDim == 3
        out(nn).Z = zeros(1,numel(TrackLinksnn));
        out(nn).std_z = zeros(1,numel(TrackLinksnn));
        out(nn).sigma_z = zeros(1,numel(TrackLinksnn));
        out(nn).std_sigma_z = zeros(1,numel(TrackLinksnn));
    end
    if wvDim
        out(nn).Wv = zeros(1,numel(TrackLinksnn));
        out(nn).std_wv = zeros(1,numel(TrackLinksnn));
        out(nn).sigma_wv = zeros(1,numel(TrackLinksnn));
        out(nn).std_sigma_wv = zeros(1,numel(TrackLinksnn));
    end
    if isfield(movieInfo,'roiIdxAll')
        out(nn).roiIdxAll = zeros(1,numel(TrackLinksnn));
    end
    if isfield(movieInfo,'roiIdxFrame')
        out(nn).roiIdxFrame = zeros(1,numel(TrackLinksnn));
    end
    if isfield(movieInfo,'model')
        out(nn).model = zeros(1,numel(TrackLinksnn));
    end
    if isfield(movieInfo,'modelIdx')
        out(nn).modelIdx = zeros(1,numel(TrackLinksnn));
    end
    for tt = TrackLinksnn
        count = count+1;
        out(nn).pH1(count) = movieInfo(tt).pH1(stats.TrackLinks(nn,tt),1);
        %get x info
        out(nn).X(count) = movieInfo(tt).xCoord(stats.TrackLinks(nn,tt),1);
        out(nn).sigma_x(count) = movieInfo(tt).xSigma(stats.TrackLinks(nn,tt),1);
        %get y info
        out(nn).Y(count) = movieInfo(tt).yCoord(stats.TrackLinks(nn,tt),1);
        out(nn).sigma_y(count) = movieInfo(tt).ySigma(stats.TrackLinks(nn,tt),1);
        %get z info
        if probDim == 3
            out(nn).Z(count) = movieInfo(tt).zCoord(stats.TrackLinks(nn,tt),1);
            out(nn).sigma_z(count) = movieInfo(tt).zSigma(stats.TrackLinks(nn,tt),1);
        end
        %get wv info
        if wvDim
            out(nn).Wv(count) = movieInfo(tt).wvCoord(stats.TrackLinks(nn,tt),1);
            out(nn).sigma_wv(count) = movieInfo(tt).wvSigma(stats.TrackLinks(nn,tt),1);
        end
        %get photons info
        out(nn).Photons(count) = movieInfo(tt).photons(stats.TrackLinks(nn,tt),1);
        %get fitting info
        if isfield(movieInfo,'roiIdxAll')% & ~isempty(movieInfo(tt).roiIdxAll)
            out(nn).roiIdxAll(count) = movieInfo(tt).roiIdxAll(stats.TrackLinks(nn,tt));
        end
        if isfield(movieInfo,'roiIdxFrame')%& ~isempty(movieInfo(tt).roiIdxFrame)
            out(nn).roiIdxFrame(count) = movieInfo(tt).roiIdxFrame(stats.TrackLinks(nn,tt));
        end
        if isfield(movieInfo,'model')
            out(nn).model(count) = movieInfo(tt).model(stats.TrackLinks(nn,tt));
        end
        if isfield(movieInfo,'modelIdx')
            out(nn).modelIdx(count) = movieInfo(tt).modelIdx(stats.TrackLinks(nn,tt));
        end
        %document frame
        if isfield(movieInfo,'frame')
            out(nn).Frame(count) = movieInfo(tt).frame+1;
        else
            out(nn).Frame(count) = tt;
        end
        if size(movieInfo(tt).xCoord,2)> 1
            %get x info
            out(nn).std_x(count) = movieInfo(tt).xCoord(stats.TrackLinks(nn,tt),2);
            out(nn).std_sigma_x(count) = movieInfo(tt).xSigma(stats.TrackLinks(nn,tt),2);
            %get y info
            out(nn).std_y(count) = movieInfo(tt).yCoord(stats.TrackLinks(nn,tt),2);
            out(nn).std_sigma_y(count) = movieInfo(tt).ySigma(stats.TrackLinks(nn,tt),2);
            %get z info
            if probDim == 3
                out(nn).std_z(count) = movieInfo(tt).zCoord(stats.TrackLinks(nn,tt),2);
                out(nn).std_sigma_z(count) = movieInfo(tt).zSigma(stats.TrackLinks(nn,tt),2);
            end
            %get wv info
            if wvDim
                out(nn).std_wv(count) = movieInfo(tt).wvCoord(stats.TrackLinks(nn,tt),2);
                out(nn).std_sigma_wv(count) = movieInfo(tt).wvSigma(stats.TrackLinks(nn,tt),2);
            end
            %get photons info
            out(nn).std_Photons(count) = movieInfo(tt).photons(stats.TrackLinks(nn,tt),2);    
        end
    end
end

if ~exist('out','var')
    error('no trajectories were created, check parameters');
end
       

stats.timeGapCloseLAP = toc(t);
fprintf('%s%g%s\n','Elapsed time is ',stats.timeGapCloseLAP,' seconds.')

[tracks, tracksStd, subRegionModel] = getTracks(out,costMatOptions.dataSize);
