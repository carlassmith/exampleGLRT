function params = getDefaultParamsGeneral
% getDefaultParamsGeneral get default values for ParamsGeneral
%
% Example:
%   params = getDefaultParamsGeneral;
%
% ParamsGeneral:  General parameters for data loading and preprocessing
%
% Fields:
%       Verbose - (binary) showing progress bars, printing to screen
%             etc. Default 1.
%       SaveAppendix - (string) string to append on the end of
%             the DataBaseName for saving. if empty just use 
%             SaveName and DataBaseName are identical 
%       Frames - (1 by nframes array) frames to analyze (1 based). if empty to analyze all
%            frames. Default []
%       Psf - (scalar) expected sigma for psf. Default 1
%       Roi - (1 by 4 array) region of interest to [xmin xmax ymin ymax]. If empty
%         use full image. Default []
%       Gain - (scalar) gain factor which scales (division) the data to be poisson
%          distributed. Default 1
%       CCDOffset - (scalar) constant background to subtract from the data. Default 0.
%       VariableName - (string) name of variable in which data is saved in
%                  DataFile. Default 'sequence'
%       PixelSize - (scalar) pixel size in microns. If empty display in microns
%               is not available. Default []
%       TimeStep - (scalar) frame rate in seconds. If empty display in seconds
%               is not available. Default []
%
% Created by Pat Cutler January 2012 (UNM)

params = struct('Verbose',1,'SaveAppendix','','Frames',[],'Psf',1,...
    'Roi',[],'Gain',1,'CCDOffset',0,'VariableName','sequence',...
    'PixelSize',[],'TimeStep',[]);


end

