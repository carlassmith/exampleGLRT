function [CRLB I] = fisherInformationGauss(Npixels,sigma,I0,bg,offset)
%
% %     |xx  xy  xI  xbg  |
% % I = |yx  yy  yI  ybg  |
% %     |Ix  Iy  II  Ibg  |
% %     |bgx bgy bgI bgbg |
%
%     I = 1;
%     sigma = 30;
%     bg = 1;

intpol = 'linear';
if((I0 == 0)||(sigma == 0)||(Npixels == 0))
    error('Arguments Npixels, sigma and I have to be bigger then 0.');
    return;
end
if (bg == 0)
    bg = 10^-10;
end

if nargin<5
    offset=[0 0];
end

im = newim([Npixels, Npixels]);

theta_x=offset(1);
theta_y=offset(2);

X=xx([Npixels Npixels]);
X=X-min(X);
X=single(X-(Npixels-1)/2);
Y=yy([Npixels Npixels]);
Y=Y-min(Y);
Y=single(Y-(Npixels-1)/2);

PSFx=1/2*(erf((X-theta_x+.5)/(sqrt(2)*sigma))-erf((X-theta_x-.5)/(sqrt(2)*sigma)));
PSFy=1/2*(erf((Y-theta_y+.5)/(sqrt(2)*sigma))-erf((Y-theta_y-.5)/(sqrt(2)*sigma)));
mu_k=I0*PSFx.*PSFy+bg;

dmudx_k=I0*sqrt(2/pi)/2/sigma*(exp(-1/2*((X+.5-theta_x)/sigma).^2)-...
    exp(-1/2*((X-.5-theta_x)/sigma).^2)).*PSFy;

dmudy_k=I0*sqrt(2/pi)/2/sigma*(exp(-1/2*((Y+.5-theta_y)/sigma).^2)-...
    exp(-1/2*((Y-.5-theta_y)/sigma).^2)).*PSFx;

dmudI_k=PSFx.*PSFy;

dmudbg_k=1;

clear I
%Ix,x
I(1,1) = sum(sum(dmudx_k.*dmudx_k./mu_k));
%Iy,y
I(2,2) = sum(sum(dmudy_k.*dmudy_k./mu_k));
%II,I
I(3,3) = sum(sum(dmudI_k.*dmudI_k./mu_k));
%Ibg,bg
I(4,4) = sum(sum(dmudbg_k.*dmudbg_k./mu_k));

%Ix,y
I(1,2) = sum(sum(dmudx_k.*dmudy_k./mu_k));
I(2,1) = I(1,2);

%Ix,I
I(1,3) = sum(sum(dmudx_k.*dmudI_k./mu_k));
I(3,1) = I(1,3);

%Ix,bg
I(1,4) = sum(sum(dmudx_k.*dmudbg_k./mu_k));
I(4,1) = I(1,4);

%Iy,I
I(2,3) = sum(sum(dmudy_k.*dmudI_k./mu_k));
I(3,2) = I(2,3);

%Iy,bg
I(2,4) = sum(sum(dmudy_k.*dmudbg_k./mu_k));
I(4,2) = I(2,4);

%II,bg
I(3,4) = sum(sum(dmudbg_k.*dmudI_k./mu_k));
I(4,3) = I(3,4);

I = ((I > 10^-10) | (I < -10^-10)).*I;
CRLB = I^(-1);
CRLB = ((CRLB > 10^-10) | (CRLB < -10^-10)).*CRLB;

