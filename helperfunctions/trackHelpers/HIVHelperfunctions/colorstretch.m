function cm = colorstretch(value,min_max,cmap)
% COLORSTRETCH    make color map to stretch across given range
% 
% cm = colorstretch(value,min_max)
% 
% INPUTS
%   value - value within range to stretch. Can be scalar or vector.
%   min_max - min and max of the range of values to stretch color map
%   cp - colormap to use for stretching colors. Needs to be an N by 3 
%        matrix with values between 0 and 1.  Default value is rgbstretch(0:100,[0 100]). 
% OUTPUT
%   cm - color map for value. Size is length(value) by 3. Columns represetn
%        red, green, and blue respectively. Values closer to min will be
%        represented as blues and values closer to max will be represented
%        as reds.
% 
% Created by Pat Cutler March 2011 (UNM)
% 
% if nargin < 2 || isempty(min_max)
%     min_max = [500 850];
% end
% if nargin < 3 || isempty(cmap)
%     cmap = rgbstretch(0:100,[0 100]);
% end
% 
% % if sum(value < min_max(1)) || sum(value > min_max(2))
% %     error('value must be within range of min_max')
% % end
% 
% if min_max(1) == min_max(2)
%     cm = cmap(round(length(cmap)/2),:);
% else
%     value0 = round((value-min_max(1))/diff(min_max)*(size(cmap,1)-1))+1;
%     value0(value0 < 1) = 1;
%     value0(value0 > size(cmap,1)) = size(cmap,1);
%     cm = cmap(value0,:);
% end
% 
cm = cmap(value,:);
