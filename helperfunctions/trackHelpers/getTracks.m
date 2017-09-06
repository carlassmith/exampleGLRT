function [tracks tracksStd subRegionModel] = getTracks(out,dataSize)
% GETTRACKS     get reformatted trajectories and localization accuracies
%
% DESCRIPTION:
%   Returns track position in the format numTracks by time by dimensions
%
% EXAMPLE:
%   [tracks tracksStd subRegionModel] = getTracks;
%
% OUTPUTS:
%   tracks - (array) trajectories
%   tracksStd - (array) theoretical localization error for trajectories
%   subRegionModel - (array) dim 3 is roiIdxAll, roiIdxFrame, model, modelIdx
%
% Created by Pat Cutler January 2012 (UNM)

if isempty(out)
    obj.connectFits;
end


if isfield(out,'Z')
    ndims = 3;
else
    ndims = 2;
end


tracks = zeros(numel(out),dataSize(3),ndims);
tracksStd = zeros(numel(out),dataSize(3),ndims);
subRegionModel = zeros(numel(out),dataSize(3),4);


for nn = 1:numel(out)
    tracks(nn,out(nn).Frame,1:2) = [out(nn).X' out(nn).Y'];
    tracksStd(nn,out(nn).Frame,1:2) = [out(nn).std_x' out(nn).std_y'];
    if isfield(out(nn),'roiIdxAll')
        subRegionModel(nn,out(nn).Frame,:) = [out(nn).roiIdxAll'...
            out(nn).roiIdxFrame' out(nn).model' out(nn).modelIdx'];
    end
end