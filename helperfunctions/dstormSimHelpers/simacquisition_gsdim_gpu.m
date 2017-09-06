% Simacquisition_gsdim_gpu  Simulate a 2d localization microscopy acquisition using the GPU
%
% function [image_out, emitter_states] = simacquisition_gsdim_gpu(object,sample_params,optics_params,camera_params,acquisition_params)
%
function varargout = simacquisition_gsdim_gpu(varargin)

d = struct('menu','FRC resolution',...
           'display','GSDIM simulation',...
           'inparams',struct('name',       {'object_image',      'sample_params',     'optics_params',      'camera_params',     'acquisition_params'},...
                             'description',{'Object image',      'Sample parameters', 'Optics parameters',  'Camera parameters', 'Acquisition parameters'},...
                             'type',       {'image',             'cellarray',         'cellarray',          'cellarray',         'cellarray'},...
                             'dim_check',  {2,                   [],                  [],                   [],                  []},...
                             'range_check',{[],                  [],                  [],                   [],                  []},...
                             'required',   {1,                   0,                   0,                    0,                   0},...
                             'default',    {[],                  {[]},                {[]},                 {[]},                {[]}}...
                              ),...
           'outparams',struct('name',{'image_out',                      'emitter_states'},...
                              'description',{'Simulated acquisition',   'True emitter states'},...
                              'type',{'image',                          'array'}...
                              )...
           );       
       
if nargin == 1
   s = varargin{1};
   if ischar(s) && strcmp(s,'DIP_GetParamList')
      varargout{1} = struct('menu','none');
      return
   end
end

try
    [object_im, ~, ~, ~, ~] = getparams(d,varargin{1});
    sample_params = varargin{2};
    optics_params = varargin{3};
    camera_params = varargin{4};
    acquisition_params = varargin{5};    
catch
    if ~isempty(paramerror)
        error(paramerror)
    else
        error('Parameter parsing was unsuccessful.')

    end
end

if any(object_im<0) || any(isnan(object_im)) || any(isinf(object_im))
    error('Illegal values in object image.');
end

%% Input parameters

% Sample properties
N_emitters = sample_params.N;                   % Expected number of emitters
bg = 0;                          % Background [photons/pixel]                 

% Fluorophore properties
% Emitter switching rate constants
N_emitters = sample_params.N;                   % Expected number of emitters
k_on = sample_params.k_on;
k_off = sample_params.k_off;
k_b = sample_params.k_b;

n_photons = 1;              % Expected collected photons per full frame
lambda = sample_params.lambda;                  % Emission wavelength [nm]

% Properties of optical system
NA = optics_params.NA;                          % Numerical aperture;

% Camera properties
sz = camera_params.fov;                         % Field of view [CCD pixels]
CCD_pixelsize = camera_params.ps;               % Backprojected camera pixel size [nm]
CCD_offset = camera_params.offset;              % Pixel offset
gain = camera_params.gain;                      % Gain: ADU / photon
readnoise = camera_params.readnoise;            % Readout noise [ADU/pixel]
CCD_precision = camera_params.precision;        % Data format precision

% Acquisition parameters
timeframes = acquisition_params.timeframes;     % Duration of the acquisition
% subframes = acquisition_params.subframes;       % Subframes for emitter switching

%% Compute parameters
PSFsigma = 0.3*lambda/NA/CCD_pixelsize;         % Standard deviation of the Gaussian PSF

%% Draw emitter positions
object_im(object_im<0) = 0;
object_im = object_im./sum(object_im)*N_emitters;
object_im = noise(object_im,'poisson');
emitter_locs = gotopoints(object_im);

% Rescale positions to match the CCD field of view
emitter_locs(:,1) = emitter_locs(:,1)*(sz(1)/imsize(object_im,1));
emitter_locs(:,2) = emitter_locs(:,2)*(sz(2)/imsize(object_im,2));
N_emitters = size(emitter_locs,1);

% %% Switching simulation
% 
% % Define Markov model
% trans_mat = [(1-k_on/subframes) k_on/subframes 0; k_off/subframes 1-(k_off+k_b)/subframes k_b/subframes; k_b 0 1];
% visible_mat = [1 0;0 1;1 0];
% 
% % Simulate switching
% emitter_states = [];
% 
% for nn = 1:N_emitters
%     [seq states] = hmmgenerate(timeframes*subframes,trans_mat,visible_mat);
%     seq = seq - 1;
%     seq = sum(reshape(seq,[subframes timeframes]),1);
%     % Format of emitter_states: [emitter_id x y t n_photons]
%     if (sum(seq~=0) ~= 0)
%         emitter_states = cat(1,emitter_states,[ones(sum(seq~=0),1)*[nn emitter_locs(nn,:)] find(seq)' nonzeros(seq)]);
%     end
% end
% 
% emitter_states(:,5) = emitter_states(:,5)*n_photons/subframes;
% emitter_states = sortrows(emitter_states,4);
% 
% clear seq states

%% Model based state simulation

emitter_states = [];    % Format of emitter_states: [emitter_id x y t n_photons]

% Define Markov model
k_act = (k_on)*(k_off+k_b)/(k_on+k_off+k_b);
k_bl = k_act - k_on*k_off/(k_on+k_off+k_b);
k_off = (k_off+k_b);

% Simulate switching

M_poisson = poissrnd(k_act*timeframes,N_emitters,1);
M_geo = 1+geornd(k_bl/k_act,N_emitters,1);
M = min(M_poisson,M_geo);


for nn = 1:N_emitters
    % Activation frame numbers of emitter nn
    t = sort(randperm(timeframes,M_poisson(nn))); 
    t = t(1:M(nn));
        
    t_start = rand(M(nn),1);        	% Fractional starting time of activations
    dt = exprnd(1/(k_off+k_b),M(nn),1); % Durations of on-events
    M_subs = ceil(t_start+dt);           % Number of frames per on-event
    M_max = max(M_subs);
    
    % For every activation event
    for ll = 1:M(nn)
        if M_subs(ll)==1
            I_vec = dt(ll);
            t_vec = t(ll);
        else
            % Vector of intensities during each frame of activation event ll
            I_vec = ones(M_subs(ll),1);
            I_vec(1) = 1-t_start(ll);
            I_vec(end) = rem(t_start(ll)+dt(ll),1);
            % Vector of frame numbers of frames in activation event
            t_vec = t(ll)-1+(1:M_subs(ll))';
        end
        emitter_states = cat(1,emitter_states,[ones(M_subs(ll),1)*[nn emitter_locs(nn,:)] t_vec I_vec]);
    end
end

emitter_states(:,5) = emitter_states(:,5)*n_photons;
emitter_states = sortrows(emitter_states,4);
emitter_states = emitter_states(emitter_states(:,4)<=timeframes,:);

%% Create movie
im_out = newim(sz(1),sz(2),timeframes);

% Find start and end indices of entries in emitter_states for each time
% point
ed = 1:size(emitter_states,1)-1;
ed = ed(emitter_states(2:end,4)>emitter_states(1:end-1,4));
ed = ed';
st = [1; ed+1];
ed = [ed; size(emitter_states,1)];

% Generate noisefree images on GPU
for ii = 1:length(st)
   sz = single(sz);
   xy_tmp = single(emitter_states(st(ii):ed(ii),2:3));
   t_tmp = emitter_states(st(ii),4);
   I_tmp = single(emitter_states(st(ii):ed(ii),5));
   s_tmp = PSFsigma*I_tmp./I_tmp;
   
   im_tmp = GPUgenerateBlobs(sz(1),sz(2),xy_tmp(:,1),xy_tmp(:,2),I_tmp,s_tmp,s_tmp,single(zeros(size(xy_tmp,1),1)),1);
   im_out(:,:,emitter_states(st(ii),4)-1) = im_tmp;
end

varargout{1}=im_out;
varargout{2} = emitter_states;