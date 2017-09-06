function params = getDefaultParamsFit
% getDefaultParamsFit get default values for ParamsFit
%
% Example:
%   params = getDefaultParamsFit;
%
% ParamsFit: Localization Parameters 
%
% Fields:
%       FitType - (integer) specifies FitType for GPUgaussMLEv2.
%           If empty fit with FitType 1 to estimate position 
%           followed by FitType 2 to etimate PSF sigma. Default [].
%       BoxSize - (scalar) fit box size in pixels. Default 7
%       Iterations - (scalar) iterations of Newton's method for 
%           GPUgauss fitting. Default 10
%       MaxCudaFits - (scalar) max cuda fits. Default 100000
%

 params = struct('FitType',[],'Iterations',10,'MaxCudaFits',100000,...
     'BoxSize',7);
             

end

