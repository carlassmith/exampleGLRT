function [gausint bg] = finiteGaussPSFerf(Npixels,sigma,I,bg,cor,aX,aY)
%%
% [gausint bg] = finiteGaussPSFerf(Npixels,sigma,I,bg,cor,aX,aY)
    if (nargin < 5)
        cor = [(Npixels-1)/2 (Npixels-1)/2 0];
    end
    if (nargin < 7)
        aX=0;
        aY=0;
    end

%     if((I == 0)||(sigma == 0)||(Npixels == 0))
%         error('Arguments Npixels, sigma and I have to be bigger then 0.');
%         return;
%     end
    if (bg == 0)
        bg = 10^-10;
    end    
    
    Ncor=size(cor,1);
    
    X=yy([Npixels Npixels Ncor]);
    X=single(X-min(X));
    Y=xx([Npixels Npixels Ncor]);
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
    
    ii = repmat(xx([Npixels Npixels]),[1 1 size(Xpos,3)]);
    ii=ii+abs(min(ii));
    jj = repmat(yy([Npixels Npixels]),[1 1 size(Xpos,3)]);
    jj=jj+abs(min(jj));
    if min(bg+aX.*ii+aY.*jj) < 0
       bg=bg+abs(min(bg+aX.*ii+aY.*jj));
       disp('Warning bg changed!');
    end
    
    
    gausint=I/4.*((erf((X-Xpos+.5)./(sqrt(2)*sigmax))-erf((X-Xpos-.5)./(sqrt(2)*sigmax))).*...
        (erf((Y-Ypos+.5)./(sqrt(2)*sigmay))-erf((Y-Ypos-.5)./(sqrt(2)*sigmay))))+bg+aX.*ii+aY.*jj;

%     gausint=dip_image(gausint);
end    
    
    
    
    
    
