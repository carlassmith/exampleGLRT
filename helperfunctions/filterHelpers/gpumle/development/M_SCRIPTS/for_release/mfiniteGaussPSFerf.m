%finiteGaussPSFerf   Make Gaussian Spots using finite pixel size
%
%   [out] = mfiniteGaussPSFerf(Npixels,sigma,I,bg,cor,avg,flag,max)
%
%       INPUT
%   Npixels:    linear size in pixels
%   sigma:      PSF sigma in pixels, scaler gives symmetric, [sx sy] gives
%               asymmetric.
%   I:          Photons/frame
%   bg:         background/pixel
%   cor:        coordinates array size of [2 N] where N is number of spots 
%               to generate.
%   avg:        average number of spots in an image, default = 1
%   flag:       if set to 0, spots are randomly drawn about the avg. if set
%               to 1, there are always an avg. number of spots, default = 1
%   mxe:        maximum number of emitters in a given frame, default = 8
%
%       OUTPUT
%   out:        3D stack of images. 


function [out] = mfiniteGaussPSFerf(Npixels,sigma,I,bg,cor)%,avg,flag,mxe)


calcInt=1;

% if (nargin < 8)
%     mxe = 8;
% end
% if (nargin < 7)
%     flag = 1;
% end
% if (nargin < 6)
%     avg = 1;
% end
    
if (nargin < 5)
    cor = [(Npixels-1)/2 (Npixels-1)/2 0];
end
if (nargin <4)
   error('Minimal usage: finiteGaussPSFerf(Npixels,sigma,I,bg)');
end

if (bg == 0)
    bg = 10^-10;
end    

Ncor=size(cor,1);

X=xx([Npixels Npixels Ncor]);
X=single(X-min(X));
Y=yy([Npixels Npixels Ncor]);
Y=single(Y-min(Y));

Xpos=repmat(shiftdim(cor(:,1),-2),[Npixels,Npixels,1]);
Ypos=repmat(shiftdim(cor(:,2),-2),[Npixels,Npixels,1]);

if size(I,1) > 1
    I=repmat(shiftdim(I,-2),[Npixels,Npixels,1]);
    if max(size(Xpos) ~= size(I))
         error('Size of I and maybe others are incorrect.');
    end
end
if size(bg,1) > 1
    bg=repmat(shiftdim(bg,-2),[Npixels,Npixels,1]);
    if max(size(Xpos) ~= size(bg))
         error('Size of bg and maybe others are incorrect.');
    end
end

if length(sigma)==1
    sigmay = sigma;
    sigmax = sigma;
else    
    sigmax = sigma(1);
    sigmay = sigma(2);
end

if calcInt
gausint=I/4.*((erf((X-Xpos+.5)./(sqrt(2)*sigmax))-erf((X-Xpos-.5)./(sqrt(2)*sigmax))).*...
    (erf((Y-Ypos+.5)./(sqrt(2)*sigmay))-erf((Y-Ypos-.5)./(sqrt(2)*sigmay))))+bg;
else
   gausint = I/(2*pi*sigmax*sigmay)*exp(-1/2*((X-Xpos)/sigmax).^2).*exp(-1/2*((Y-Ypos)/sigmay).^2)+bg;    
end

% add multiple gaussians together depending on flag type, avg, and max
% mod by PKR, UNM December 2012

out = gausint;

%out=dip_image(gausint);


end
