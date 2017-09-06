% GOTOPOINTS Go from a binned image to point locations
%     POSMAT = gotopoints(IM) converts the binned localization data  in
%     IM to point location data in the N-by-2 array posmat by
%     taking the number of points in each pixel bin and associating them with random
%     locations in the pixel area.

function out = gotopoints(varargin)

d = struct('menu','Localization Microscopy',...
           'display','Go to points',...
           'inparams',struct('name',       {'in'},...
                             'description',{'Image with binned positions'},...
                             'type',       {'image'},...
                             'dim_check',  {0},...
                             'range_check',{[]},...
                             'required',   {1},...
                             'default',    {[]}...
                              ),...
           'outparams',struct('name',{'out'},...
                              'description',{'List of positions'},...
                              'type',{'array'}...
                              )...
           );       


if nargin == 1
   s = varargin{1};
   if ischar(s) & strcmp(s,'DIP_GetParamList')
      out = d;
      return
   end
end

try
   in = getparams(d,varargin{:});
catch
   if ~isempty(paramerror)
      error(paramerror)
   else
      error('Parameter parsing was unsuccessful.')
   end
end

if sum(abs(round(in)-in))~=0
    error('Input image can only contain integer values.')
end

% 
in = im2mat(in);
out = [];

for nn = 1:max(max(in))
    % Find positions of nonzero pixels and add locations to out
    [y x] = find(in);
    out = cat(1,out,[x y]);
    % Subtract 1 from each nonzero pixel
    in = in - (in>0);
end

% Add random displacement of position within the pixel
out = out - 1.5*ones(size(out))+rand(size(out));