function params = getDefaultParamsPreFilterFits(fitSigma)
% getDefaultParamsFilterFits get default values for ParamsFilterFits
%
% Example:
%   params = getDefaultParamsFilterFits;
%
% ParamsFilterFits: Localization Parameters
%
% Fields:
%   circularityMax
%   circularityMin
%   pH1Max
%   pH1Min

params.circularityMax = 3;
params.circularityMin = 0.5;
params.pH1Max = 1;
params.pH1Min = 0;
params.minPixelDist = 0;
params.clusterSizeMin = 0;
params.clusterSizeMax = 100;

if nargin > 0 && fitSigma
    paramsPreFilterFits.MaxSigmax = Inf;
    paramsPreFilterFits.MinSigmax = 0;
    paramsPreFilterFits.MaxSigmay = Inf;
    paramsPreFilterFits.MinSigmay = 0;
end

end