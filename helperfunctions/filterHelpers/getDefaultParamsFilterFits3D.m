function params = getDefaultParamsFilterFits
% getDefaultParamsFilterFits get default values for ParamsFilterFits
%
% Example:
%   params = getDefaultParamsFilterFits;
%
% ParamsFilterFits: Localization Parameters
%
% Fields:
%       MaxCRLBSTD - (scalar) Maximum CRLB for position parameters (x,y). Default 0.5
%       MinPValue -  (scalar) minimum probability value from gamma
%                    distribution. Default 0.01
%       MinPhotons - (scalar) minimum non-zero photon. Default 50. 
%       MinPixelDist - (scalar) minimum pixel distance allowed
%           between two localized objects. Default 3.

params = struct('MaxCRLBSTD',0.5,'MinPValue',0.01,...
    'MinPhotons',50,'MinPixelDist',3);


end

