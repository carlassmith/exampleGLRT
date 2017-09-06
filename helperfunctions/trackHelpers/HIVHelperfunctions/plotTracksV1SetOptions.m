function options = plotTracksV1SetOptions(varargin)
% plotTracksV1SetOptions   set default options for plotTracksV1
% 
% options = plotTracksV1SetOptions;
%
% OUTPUT:
%   options: structured array with the following fields:
%
%       tracks: trajectories with dims [tracks by time 2]. Default []
%       h: figure handle. Default []
%       ha: axes handle. Default []
%       colorByValue: binary for coloring tracks by tracksVal. Default 1.
%       tracksVal: track value for coloration and data cursor. Default []
%       cmap: color map to use for coloring tracks. Default jet(100)
%       colorbar: binary for including colorbar. Default 1.
%       xlim: limits for x axis. Default []
%       ylim: limits for y axis. Default []
%       t: frame rate or time between consecutive observations. Default 1
%       linewith: width of lines. Default 1.5
%       linestyle: style of lines. Default '-'
%       valueName: name used for colorbar title. Default '\lambda value'
%       trackNum: track number(s) to plot. Default []
%       view: view angle for 3D axes see view. Default [45 45]
%       color: color to use for tracks if colorByValue is 0. Default [0 0 0]
%       minMaxVal: min and max value for coloring tracks.
%       markersize: size of marker identifying vertices. Default 3
%       marker: shape of marker for each vertex. Default 'o'
%       trange: time range to plot in frames. If empty plot all frames.
%           Default []
%       addShadow: relative width of shadow in comparison with linewidth
%           field. Default 0.5
%       highlightTracks: binary 1- highlight tracks when selected with data
%           cursor. 0 - don't highlight tracks. Default 1
%       displayName: string describing patch object. Default ''
%
%       
% NOTE: define vert, face, and val options at your own risk. The
% respective dimensions must correspond for building the patch object.
% 
% See also plotTracksV1
% 
% Created by Pat Cutler January 2010 (UNM)

options = struct;

if (nargin == 0 || ischar(varargin{1}))
    options = struct('tracks',[],...
        'h',[],...
        'ha',[],...
        'colorByValue',1,...
        'tracksVal',[],...
        'cmap',jet(100),...
        'colorbar',1,...
        'xlim',[],...
        'ylim',[],...
        't',1,...
        'linewidth',1.5,...
        'linestyle','-',...
        'valueName','\lambda value',...
        'trackNum',[],...
        'view',[45 45],...
        'color',[0 0 0],...
        'minMaxVal',[],...
        'markersize',3,...
        'marker','o',...
        'trange',[],...
        'addShadow',.5,...
        'highlightTracks',1,...
        'displayName','');
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