%% Simulation input parameters
% object_im =  cos(pi*rr(128,128)/128)*cos(phiphi(128,128)*10).^2;
% object_im = (object_im>0)*object_im;
% 
object_im = newim(128,128);
object_im(:,:) = 1.0;

sample_params = [];
optics_params = [];
camera_params = [];
acquisition_params = [];  
    
% Acquisition parameters
acquisition_params.timeframes = 5E4;    % Duration of the acquisition
acquisition_params.subframes = 10;      % Subframes for emitter switching
acquisition_params.fitbox = 9;          % Size of ROIs for single emitter fitting [CCD pixels]

% Sample properties
sample_params.N_sites = 2000;          % Number of binding sites
sample_params.N_persite = 0;            % Number of emitters per binding site
sample_params.labelsize = 0.05;         % Size of fluorescent labels [nm]
sample_params.bg = 1;                   % Background [photons/pixel]                 

% Fluorophore properties
% Emitter switching rate constants
sample_params.k_on = 5E-4;
sample_params.k_off = 0.5;
sample_params.k_boff = 0.1*sample_params.k_on;
sample_params.k_bon = 0.01*sample_params.k_off;


% Emission wavelength
sample_params.photons = 1000;           % Expected collected photons per full frame
sample_params.lambda = 670;             % Wavelength [nm]

% Properties of optical system
optics_params.NA = 1.45;
camera_params.fov = imsize(object_im);  % Field of view [CCD pixels]

% Camera properties
camera_params.ps = 100;                 % Pixel size [nm]

%% Run simulation
tic;
[coords, emitter_states] = simgsdim_fitdata_model(object_im,sample_params,optics_params,camera_params,acquisition_params);
toc

%% Q stats
M = hist(emitter_states(:,1),1:max(emitter_states(:,1)));
Q_gt = mean(M.*(M-1))./mean(M);