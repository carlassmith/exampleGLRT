function [ im_out, emitter_states, GTImage, numOfTrueP ] = createDStormSim( options )
%createDStormSim Creates DStorm Simlation
% SYNOPSIS:
%   [ image_out, emitter_states, GTImage, numOfTrueP,imgOutNoiseFree] = createDStormSim( options )
% PARAMETERS:
%    options = dipDStormGetDefaultParams(imageSize)
% 
% OUTPUTS:
%   image_out: Simlated DSTORM with noise
%   emitter_states: emitter_states (FORMAT [emitter_id x y t n_photons])
%   GTImage: pixelated locations
%   numOfTrueP: Number of spots
%   imgOutNoiseFree: imlated DSTORM without noise
% 
%   SEE ALSO:
%       dipDStormGetDefaultParams
[im_out, emitter_states] = simacquisition_gsdim_gpu(options.object_im,options.sample_params,options.optics_params,options.camera_params,options.acquisition_params);


% emitter_states = [emitter_id x y t n_photons]
XYT = emitter_states(:,2:4);
GTImage = coord2image([XYT(:,2),...
    XYT(:,1),...
    XYT(:,3)-1],...
    imsize(im_out));

%% TODO: or just check size?
numOfTrueP = sum(GTImage);
end

