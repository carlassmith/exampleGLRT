function options = dipDStormGetDefaultParams
%dipDStormGetDefaultParams Gets Defaults DStorm Simulation parameters
% SYNOPSIS:
%   options = dipDStormGetDefaultParams
% PARAMETERS:
%    -
% OUTPUTS:
%   options: default settings struct
% 
% TODO: make it into a setfield function!! Now only returns default params

%% Input parameters
object_im = [];
sample_params = [];
optics_params = [];
camera_params = [];
acquisition_params = [];  
    
% Acquisition parameters
acquisition_params.timeframes = 1000;                   % Duration of the acquisition

% Sample properties
sample_params.N = [];               % Number of emitters
sample_params.bg = [];                % Background [photons/pixel]                 

% Fluorophore properties
% Emitter switching rate constants
sample_params.k_on = 0.1;
sample_params.k_off = 1;
sample_params.k_b = 0.1*sample_params.k_off;

% Emission wavelength
sample_params.photons = 1000;       % Expected collected photons per full frame
sample_params.lambda = 670;         % Wavelength [nm]

% Properties of optical system
optics_params.NA = 1.45;


% Camera properties
camera_params.ps = 100;                 % Pixel size [nm]
camera_params.offset = 0;             % Camera offset (bgin for reconstruction)
camera_params.gain = 1;
camera_params.readnoise = 0;
camera_params.precision = 'int16';

options.camera_params = camera_params;
options.sample_params = sample_params;
options.object_im = object_im;
options.optics_params = optics_params;
options.acquisition_params=acquisition_params;
end