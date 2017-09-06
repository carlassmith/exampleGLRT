function options = costMatCloseGapsSetOptions(varargin)
% COSTMATCLOSEGAPSSETOPTIONS    set gap closing connection parameters for costMatCloseGaps
% 
% options = costMatCloseGapsSetOptions;
% Creates a structure of 'options' for costMatCloseGaps set to default values.
%
% options = costMatCloseGapsSetOptions('param1',value1,'param2',value1,...);
% Creates a structure of 'options' for costMatCloseGaps setting named parameters to
% specified values and all other parameters to their default.
% 
% options = costMatCloseGapsSetOptions(oldOptions,'param1',value1,'param2',value2,...)
% Creates a copy of 'oldOptions' with the named parameters altered.
% 
% Dependendies:
% 
% See also costMatCloseGaps
% 
% Created by Pat Cutler August 2011 (UNM)

% set default values for conditions 1 and 2
if (nargin == 0 || ischar(varargin{1}))
    options = struct('funcName','costMatCloseGaps',...
        'maxSearchDistPerFrame',[2 2]',...
        'maxSearchDist',[5 5]', ...
        'maxWvSearchDist',10, ...
        'D',[0.1 0.1], ...
        'kon',0.1, ...
        'wvJump',5, ...
        'timeWindow',40,...
        'costBD',0, ...
        'blur',0.69);
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


