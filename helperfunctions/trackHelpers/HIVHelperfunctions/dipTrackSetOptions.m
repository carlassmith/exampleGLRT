function options = dipTrackSetOptions(varargin)
% dipTrackSetOptions   set options for dipTrack
% 
% options = dipTrackSetOptions;
% Creates a structure of options for dipTrack set to default values.
%
% options = dipTrackSetOptions('param1',value1,'param2',value1,...)
% Creates a structure of options for dipTrack seting named parameters to
% specified values and all other parameters to their default.
% 
% options = dipTrackSetOptions(oldopts,'param1',value1,'param2',value2,...)
% Creates a copy of oldopts with the named parameters altered.
%
% PARAMETERS:
%       im - image (required)
%       plotTracksOptions - inputs for plotTracksV1. Default plotTracksV1SetOptions
%       HeadMarkersize - 1 by 2 vector indicating the marker size for track
%                        head. Default [5 3]
%       headMarker - single character or vector of characters specifying 
%                    the shape of the head marker. If single character then
%                    a single head marker shape is used. If vector of
%                    characters is the same length as the number of tracks
%                    specified in plotTracksOptions then each track is
%                    represented with a different character. Default 'o'
%       npoints - number of tail points. Default 10.
%       h - figure handle
% 
% See also dipTrack
% 
% Created by Pat Cutler January 2012 (UNM)

options = struct;

if (nargin == 0 || ischar(varargin{1}))
    options = struct('im',[],...
        'plotTracksOptions',plotTracksV1SetOptions,...
        'npoints',10,...
        'headMarkersize',[5 3],...
        'headMarker','o',...
        'h',[],...
        'updateImage',1);
    options.plotTracksOptions.addShadow = 2;
end
if(nargin > 0)
    if(ischar(varargin{1}))
        i = 1;
        while (i < length(varargin))
            options = setfield(options,varargin{i},varargin{i+1});
            i = i+2;
        end
    else if(isstruct(varargin{1}))
            if(isstruct(varargin{1}))
                options = varargin{1};
                j = 2;
                while (j < length(varargin))
                    options = setfield(options,varargin{j},varargin{j+1});
                    j = j+2;
                end
            else disp('error - in calling function');
            end
        end
    end
end


return