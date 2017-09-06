function params = getDefaultParamsConnectFits
            % getDefaultParamsConnectFits get default values for ParamsConnectFits
            %
            % Example:
            %   params = getDefaultParamsConnectFits;
            %
            % ParamsConnectFits: Parameters for connecting localizations into tracks
            %
            %   FunctionCall - (string) evaluated for connecting
            %       localizations. This can be swapped with any user
            %       defined input. Parameters need to be adjusted accordingly. 
            %       Default 'connectCoordsLAP(obj.FitResults,obj.ParamsConnectFits);'
            %
            %     NOTE: all other parameters are dependent on the Function
            %     call. The parameters specified below are specific for
            %     connectCoordsLAP. ANY USER DEFINED FUNCTION MUST HAVE TWO
            %     OUTPUTS. The first output must be in the proper
            %     format for the Tracks property. The second output can contain
            %     any additional information that the user wants to return.
            %     The second output is written to obj.Stats.Tracks.connectFits
            %
            %   costMatF2F - (struct) parameters for making frame to
            %       frame connections. Contains the following fields:
            %
            %           funcName: (string) name of function to use to build
            %               frame to frame cost matrix. Default 'costMatFrame2FrameDensity'
            %           density: (scalar) single molecule density in
            %                       particles per frame. If empty populated
            %                       in connectCoordsLAP from localzation
            %                       density.
            %           D: (1 by 2 array) diffusion coefficients in [x y]
            %               in pixels. Default [.1 .1].
            %           maxSearchDist: (1 by 2 array) maximum search 
            %               distance for connection in [x y] in pixels.
            %               Default [2 2].
            %           kon: (scalar) characteristic rate for observing a 
            %               new probe in a subsequent frame. Default 0.1
            %           koff: (scalar) characteristic rate for not observing 
            %               the same probe in a subsequent frame. Default 0.1
            %           blur: (scalar) log prob variance allowed before connections 
            %               are rejected. Default 0.69
            %
            %   costMatGC - (struct) parameters for gap closing 
            %       connections. Contains the following fields:
            %
            %           funcName: (string) name of function to use to build
            %               gap closing cost matrix. Default 'costMatCloseGapsDensity' 
            %           density: (scalar) single molecule density in
            %                       particles per frame. If empty populated
            %                       in connectCoordsLAP from localzation
            %                       density.
            %           D: (1 by 2 array) diffusion coefficients in [x y]
            %               in pixels. Default [.1 .1].
            %           timeWindow: (scalar) maximum time to search for
            %               connection. Default 40
            %           maxSearchDistPerFrame: (1 by 2 array) maximum 
            %               search distance per frame for connection in 
            %               [x y] in pixels. Default [2 2]
            %           maxSearchDist: (1 by 2 array) maximum search 
            %               distance for connection in [x y] in pixels.
            %               Default [8 8].
            %           minTrackLen: (scalar) minimum track length for gap
            %               closing. Default 2
            %           kon: (scalar) characteristic rate for (birth) 
            %               observing a new probe in the 'timeWindow' and 
            %               for (death) not observing the same probe in the 
            %               'timeWindow'. Default 0.1
            %           blur: (scalar) log prob variance allowed before connections 
            %               are rejected. Default 0.69
            %
            % see also connectFits, connectCoordsLAP,
            % costMatFrame2FrameDensity, costMatFrame2FrameSetOptions
            % costMatCloseGapsDensity, costMatFrame2FrameSetOptions
            %
            % Created by Pat Cutler January 2012 (UNM)   
                  
            
            % The actual tracking function in use, you can swap for something else
            params.FunctionCall = 'connectCoordsLAP(obj.FitResults,obj.ParamsConnectFits);';
            
            %%% set parameters for frame 2 frame connections %%%
            params.costMatF2F = costMatFrame2FrameSetOptions;
            params.costMatF2F = rmfield(params.costMatF2F,'wvJump');
            params.costMatF2F = rmfield(params.costMatF2F,'maxWvSearchDist');
            params.costMatF2F.funcName = 'costMatFrame2FrameDensity';
            params.costMatF2F.density = [];
            params.costMatF2F.D = [.1 .1];
            params.costMatF2F.maxSearchDist = [2 2];
            params.costMatF2F.kon = 0.1;
            params.costMatF2F.koff = 0.1;
            params.costMatF2F.blur = 0.69;
            
            %%% set parameters for gap closing %%%
            params.costMatGC = costMatCloseGapsSetOptions; 
            params.costMatGC = rmfield(params.costMatGC,'wvJump');
            params.costMatGC = rmfield(params.costMatGC,'maxWvSearchDist');
            params.costMatGC.funcName = 'costMatCloseGapsDensity'; 
            params.costMatGC.density = [];
            params.costMatGC.D = [.1 .1];
            params.costMatGC.timeWindow = 40; 
            params.costMatGC.maxSearchDistPerFrame = [2 2];
            params.costMatGC.maxSearchDist = [8 8];
            params.costMatGC.minTrackLen = 2;
            params.costMatGC.kon = 0.1;
            params.costMatGC.koff = 0.1;            
            params.costMatGC.blur = 0.69;


end

