function params = getDefaultParamsFindBoxCenters
% getDefaultParamsFindBoxCenters get default values for ParamsFindBoxCenters
%
% Example:
%   params = getDefaultParamsFindBoxCenters;
%
% ParamsFindBoxCenters: Parameters for subregions of interest for localization
% 
% Fields:
%       MinPhotons - (scalar) approximate min photons for filtering.
%                    Default 40
%       unifFilter - (1 by 2 array) uniform filtering parameters. Default [2 2]
%       maxFilter - (1 by 2 array) max filtering parameters. Default [7 7]
%       fracMax - (scalar) percent of max filter to use for binary.
%                 Default []
%       filterShape - (string) filter shape for unif and maxf filters.
%                     Default 'rectangular'
%
% see also sptFindCandidates
%
% Created by Pat Cutler January 2012 (UNM)

     params = struct('MinPhotons',40); 

end

