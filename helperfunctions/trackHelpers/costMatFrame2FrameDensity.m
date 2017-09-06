function costMat = costMatFrame2FrameDensity(movieInfo,costMatParam)
% COSTMATFRAME2FRAMEDENSITY     cost matrix for linking frames
% 
% [costMat dm birth death] = costMatFrame2FrameDensity(movieInfo,costMatParam)
% 
% INPUTS
%   movieInfo - loocalization information. see connectCoordsLAP for details
%   costMatParam - structure with input parameters. Use
%                  connectCoordsCostMatSetGapMergeSplitParams to set
%                  parameters important parameters are:
% 
%                   maxSearchDist - vector with size 1 by probDim
%                                   describing the max search distance in 
%                                   each spatial dimension 
%                   D - vector with size 1 by probDim containing expected 
%                       diffusion in each spatial dimension
%                   kon - rate for blinking on
%                   koff - rate for blinking off
%                   probDim - number of spatial dimensions
%                   density - particle density. particles/pixels^2/lambda
% OUTPUTS
%   costMat - cost matrix for frame to frame linking
% 
% See Also: connectCoordsLAP, costMatFrame2FrameSetOptions
% 
% Created by Pat Cutler August 2011 (UNM)
% Modified by Carlas Smith June 2014 (UMASS)
% 
% Adapted from Jaqaman code

% function P should supply the transition probabilities as:
% P(1,1)  s_{k-1}=0 s_k=0
% P(1,2)  s_{k-1}=0 s_k=x_{k}
% P(2,1)  s_{k-1}=x_{k-1} s_k=0
% P(2,2)  s_{k-1}=x_{k-1} s_k=x_{k}
% function call
% P = transP(row,col,movieInfo,costMatParam)

%% Input

%check whether correct number of input arguments was used
if nargin ~= nargin('costMatFrame2FrameDensity')
    error('costMatFrame2FrameDensity:narginIncorrect','costMatFrame2FrameDensity: Incorrect number of input arguments!')
end

%% compute cost matrix
transPParams = struct('dist',[],'probDim',[],'LA1',[],'LA2',[],...
                      'Coords1',[],'Coords2',[]);

if isfield(movieInfo,'zCoord')
    probDim = 3;
else
    probDim = 2;
end


if ~isequal(probDim,length(costMatParam.D),length(costMatParam.maxSearchDist))
    error('costMatFrame2FrameDensity:inproperLengthDmaxSearchDistmaxSearchDistPerFrame',...
        ['costMatFrame2FrameDensity: one of the input options ''D'' or ''maxSearchDist''' ...
        ' doesn''t have the proper length of ' num2str(probDim)]);
end

fieldNames = {'xCoord' 'yCoord' 'zCoord'};
for ii = 1:probDim
dist(:,:,ii) = (createDistanceMatrix(double(eval(['movieInfo(1).' fieldNames{ii} '(:,1)'])),...
    double(eval(['movieInfo(2).' fieldNames{ii} '(:,1)'])))).^2;
end

%Assign NaN to dist corresponding to distance > searchRadius
maxSearchDist = repmat(reshape(squeeze(costMatParam.maxSearchDist)',[1 1 probDim]),[size(movieInfo(1).xCoord,1) size(movieInfo(2).xCoord,1) 1]);
dist(dist>maxSearchDist.^2) = NaN;

LA1 = zeros(size(movieInfo(1).xCoord,1) ,probDim);
LA2 = zeros(size(movieInfo(2).xCoord,1) ,probDim);

for ii = 1:probDim
    if size(movieInfo(1).xCoord,2) > 1
        LA1(:,ii) = eval(['movieInfo(1).' fieldNames{ii} '(:,2)']);
        LA2(:,ii) = eval(['movieInfo(2).' fieldNames{ii} '(:,2)']);
    end
end

kon = costMatParam.kon;
koff = costMatParam.koff;
D = costMatParam.D;
density = costMatParam.density;
distRep  = reshape(dist,[size(dist,1)*size(dist,2) size(dist,3)]);
LA1Rep =  repmat(LA1,[size(LA2,1) 1]);
LA2Rep = repmat(LA2,[size(LA1,1) 1]);

[ dm ] = transP(2,2,kon,koff,D,[],density,distRep,probDim,LA1Rep,LA2Rep,size(LA1,1),size(LA2,1));
[ birthBlock ] = transP(1,2,kon,koff,D,[],density,distRep,probDim,LA1,LA2);
[ deathBlock ] = transP(2,1,kon,koff,D,[],density,distRep,probDim,LA1,LA2);
[ lrBlock ] = transP(1,1,kon,koff,D,[],density,distRep,probDim,LA1Rep,LA2Rep,size(LA1,1),size(LA2,1));

dm = dm - real(log(movieInfo(1).pH1*movieInfo(2).pH1'));
% birthBlock=birthBlock- real(log(diag(movieInfo(2).pH1)));
% deathBlock=deathBlock- real(log(diag(1-movieInfo(1).pH1)));
lrBlock = lrBlock-real(log((1-movieInfo(2).pH1)*(1-movieInfo(1).pH1')));

%append cost matrix
costMat = [dm deathBlock;....
    birthBlock lrBlock];
costMat(isnan(costMat)) = 0;
costMat(isinf(costMat)) = 0;

%make into sparse matrix (nonlinkmarker needs to be zero)
costMat = sparse(double(costMat));

%% Running a cost matrix check to look for questionable connections
% Any connection that doesn't hold truth after a slight perturbation is
% disregarded as too uncertain to assure a good connection

[links12, ~] = lap(real(costMat),0,0,0);
mask = ones(size(dm));

tempcMat = costMat;
for ii = 1:length(links12)
    if ii <= size(mask,1)
        if links12(ii) <= size(mask,2)
            tempcMat(ii,links12(ii)) = tempcMat(ii,links12(ii)) +  costMatParam.blur;
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
    
    %append cost matrix
    costMat = [dm.*mask deathBlock; birthBlock lrBlock.*mask'];
    costMat(isnan(costMat)) = 0;
    costMat(isinf(costMat)) = 0;
    %make into sparse matrix (nonlinkmarker needs to be zero)
    costMat = sparse(double(costMat));
end
