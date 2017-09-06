function [im_out, emitter_states] = addNoiseDStormSim(im_out,emitter_states,options)

    CCD_offset = options.camera_params.offset;              % Pixel offset
    gain = options.camera_params.gain;                      % Gain: ADU / photon
    readnoise = options.camera_params.readnoise;            % Readout noise [ADU/pixel]
    CCD_precision = options.camera_params.precision;        % Data format precision

    n_photons = options.sample_params.photons;              % Expected collected photons per full frame
    %% Add noise
    emitter_states(:,5) = emitter_states(:,5)*n_photons;
    imOutNoiseFree = im_out*n_photons+ options.sample_params.bg;
    im_out = noise(imOutNoiseFree,'poisson');
    im_out = gain*im_out;
    im_out = im_out + CCD_offset;
    if ~( readnoise == 0)
        im_out = noise(im_out,'gaussian',readnoise);
    end
    im_out = dip_image(im_out,CCD_precision);

end