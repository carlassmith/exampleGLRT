% Simacquisition_gsdim_fitdata  Obtain localizations of a simulated 2d localization microscopy acquisition
%
% function [coords, emitter_states] = simgsdim_fitdata_full(object,sample_params,optics_params,camera_params,acquisition_params)
%
function varargout = simgsdim_fitdata_model(varargin)

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
N_sites = sample_params.N_sites;                % Expected number of emitters
N_persite = sample_params.N_persite;
d_label = sample_params.labelsize;              % Label size [nm]
bg = sample_params.bg;                          % Background [photons/pixel]                 

% Fluorophore properties
% Emitter switching rate constants
k_on = sample_params.k_on;
k_off = sample_params.k_off;
k_b1 = sample_params.k_boff;
k_b2 = sample_params.k_bon;

n_photons = sample_params.photons;              % Expected collected photons per full frame
lambda = sample_params.lambda;                  % Emission wavelength [nm]

% Properties of optical system
NA = optics_params.NA;                          % Numerical aperture;

% Camera properties
sz = camera_params.fov;                         % Field of view [CCD pixels]
CCD_pixelsize = camera_params.ps;               % Backprojected camera pixel size [nm]

% Acquisition parameters
timeframes = acquisition_params.timeframes;     % Duration of the acquisition
subframes = acquisition_params.subframes;       % Subframes for emitter switching
fitboxsize = acquisition_params.fitbox;         % Size of ROIs for single emitter fitting [CCD pixels]

%% Compute parameters
PSFsigma = 0.3*lambda/NA/CCD_pixelsize;         % Standard deviation of the Gaussian PSF

%% Draw emitter positions
object_im(object_im<0) = 0;
object_im = object_im./sum(object_im)*N_sites;
object_im = noise(object_im,'poisson');
emitter_locs_tmp = gotopoints(object_im);

% Rescale positions to match the CCD field of view
emitter_locs_tmp(:,1) = emitter_locs_tmp(:,1)*(sz(1)/imsize(object_im,1));
emitter_locs_tmp(:,2) = emitter_locs_tmp(:,2)*(sz(2)/imsize(object_im,2));

%% Generate multiple emitters per binding site
N_sites = size(emitter_locs_tmp,1);

emitter_locs = [];

for mm = 1:N_sites
    if N_persite > 0
        S = poissrnd(N_persite);
    else
        if N_persite < 0;
            beta = sample_params.dimerfrac;
            S = 1+binornd(1,beta);
        else
        	S = 1;
        end
    end
    
    if S>0
       emitter_locs = cat(1,emitter_locs,repmat(emitter_locs_tmp(mm,:),[S 1]));
    end
end

emitter_locs = sortrows(emitter_locs);
emitter_locs = emitter_locs + (d_label/CCD_pixelsize)*randn(size(emitter_locs));
N_emitters = size(emitter_locs,1);

%% Switching simulation

% Define Markov model
k_act = (k_on+k_b1)*(k_off+k_b2)/(k_on+k_b1+k_off+k_b2);
k_bl = k_act - k_on*k_off/(k_on+k_b1+k_off+k_b2);
k_off = (k_off+k_b2);

% Simulate switching

M_poisson = poissrnd(k_act*timeframes,N_emitters,1);
M_geo = 1+geornd(k_bl/k_act,N_emitters,1);
M = min(M_poisson,M_geo);

emitter_states = [];

for nn = 1:N_emitters
    % Format of emitter_states: [emitter_id x y t]
    if M(nn) ~= 0
        t_tmp = sort(randperm(timeframes,M_poisson(nn)))';
        emitter_states = cat(1,emitter_states,[ones(M(nn),1)*[nn emitter_locs(nn,:)] t_tmp(1:M(nn))]);
    end
end

clear seq states

%% Localization simulation

emitter_states = sortrows(emitter_states,4);

coords = zeros(size(emitter_states,1),8);
coords(:,1:3) = emitter_states(:,2:4);
coords(:,4) = 1+geornd(k_off/n_photons,size(emitter_states,1),1);
coords(:,8) = ceil(coords(:,4)/n_photons);
coords(:,6) = poissrnd(coords(:,8)*bg*fitboxsize^2)/fitboxsize^2;
coords(:,7) = PSFsigma*(1+0.02*randn(size(emitter_states,1),1));
tau = 2*pi*(coords(:,7).^2+1/12).*coords(:,6)./coords(:,4);
coords(:,5) = sqrt((coords(:,7).^2+1/12).^2./coords(:,4).*(1+4*tau+sqrt((2*tau)./(1+4*tau))));   % Localization uncertainty
coords(:,1:2) = coords(:,1:2) + (coords(:,5)*[1 1]).*randn(size(coords(:,1:2)));

coords = coords(~any(isnan(coords)'),:);

varargout{1} = coords;
varargout{2} = emitter_states;