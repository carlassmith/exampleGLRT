function params = getDefaultParamsFilterFits
% getDefaultParamsFilterFits get default values for ParamsFilterFits
%
% Example:
%   params = getDefaultParamsFilterFits;
%
% ParamsFilterFits: Localization Parameters
%
% Fields:
%       MinCRLBSTD: 0
%       MinPhotons: 0
%            MinBg: 0
%       MaxCRLBSTD: Inf
%        MaxPFAValue: 0.05
%       MaxPhotons: Inf
%            MaxBg: Inf
%     MinPixelDist: 0

params.MinCRLBSTD = 0;
% params.MinPValue = 0;
params.MinPhotons = 0;
params.MinBg = 0;

params.MaxCRLBSTD = Inf;
params.MaxPFAValue = 0.05;
params.MaxPhotons = Inf;
params.MaxBg = Inf;

params.MinPixelDist = 0;

end

