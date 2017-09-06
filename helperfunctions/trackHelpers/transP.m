


function [ pTrans ] = transP(row,col,kon,koff,D,time,density,dist,probDim,LA1,LA2,rr,cc)

% function P supplies the transition probabilities as:
% P(1,1)  s_{k-1}=0 s_k=0
% P(1,2)  s_{k-1}=0 s_k=x_{k}
% P(2,1)  s_{k-1}=x_{k-1} s_k=0
% P(2,2)  s_{k-1}=x_{k-1} s_k=x_{k}
% Created by Carlas Smith June 2014 (UMASS/TU-DELFT)

% TODO: add two sided boudndary crossing probaility for out of focus probability!!
% Karlin, S. and Taylor H. M. (1975). A First Course in Stochastic Processes, 2nd ed., Academic
% Press.
% Jacod, J. and Protter P. (2004). Probability Essentials. 2nd edition. Springer. ISBN 3-540-
% 43871-8.
% P(|W_t| > b) = 4 \infsum_k (-1)^{k} phi((2k-1)b/sqrt(t)) assuming we at
% zero.
% function p = pwt(b,t)
%     phi=@(x).5.*(1 + erf(x./sqrt(2)));
%     k = 1:100*2*D*t;
%     p = 4.*sum((-1).^(k).*phi((2.*k-1).*b./sqrt(2*D*t)));
% end

if nargin < 8
    dist = 1;
end
if nargin < 9
    probDim = 2;
end
if nargin < 7 || nargin < 7
    LA1 = 0; LA2 = 0; 
end

if nargin < 12 || isempty(rr) ||  isempty(cc)
    cc = 1; rr = 1;
end
% true = transition probabilities calculated by Pat Cutler.
classic = false;

if row == 1 & col == 1
if nargin < 4 || isempty(time)
    time = ones(size(LA2,1),1);
end
% s_{k-1}=0 s_k=0
%     lr = sum(dist,2);
    Drep = repmat(permute(D,[1 2]),[size(LA1,1) 1]); 
    timeRep = repmat(time,[1  probDim]);   
    lr = sum(dist./(2.*(2.*Drep.*timeRep+1.*(LA1.^2+LA2.^2))+eps)+0.5.*log(2.*pi.*(2.*Drep.*timeRep+1.*(LA1.^2+LA2.^2))+eps),2);
    lr = (-log(1-density.*exp(-1.*(time).*koff).*koff)).*(lr>0); % minimize costs of lr matrix
%     .*kon*exp(kon)
    pTrans = lr;   
    if isscalar(cc)
        pTrans = reshape(pTrans,[rr,cc])';
    else
        pTrans = sparse(rr,cc,pTrans)';
    end
elseif row == 1 & col == 2
% s_{k-1}=0 s_k=x_k
        if nargin < 4 || isempty(time)
            time = ones(size(LA2,1),1);
        end
        if ~classic
            Pb = ((time).*kon-log(kon));
%             +(time).*koff-log(koff)-log(density)+
        else
            Pb = (-log(density)-log(kon)-(time).*log(1-kon));
        end
        birthCost = Pb;
        birthBlock = birthCost; %lower left
        birthBlock(birthBlock==0) = NaN;
        pTrans = birthBlock;
        pTrans = diag(pTrans);
        if ~isscalar(cc)
            pTrans = sparse(pTrans);
        end        
elseif row == 2 & col == 1
% s_{k-1}=x_{k-1} s_k=0
        if nargin < 4 || isempty(time)
            time = ones(size(LA1,1),1);
        end
        
        %P(T<=delta t) = 1- exp(-koff deltat)
        if ~classic
            Pd = (eps-real(-1.*(time).*koff+log(koff)));
%             log(density)++(time).*kon-log(kon)
        else
            Pd = eps-(time).*log(1-kon);
        end
        deathCost = Pd ;
        %generate upper right and lower left block
        deathBlock = deathCost; %upper right
        deathBlock(deathBlock==0) = NaN;
        pTrans = deathBlock;
        pTrans = diag(pTrans);
        if ~isscalar(cc)
            pTrans = sparse(pTrans);
        end
elseif row == 2 & col == 2
% s_{k-1}=x_{k-1} s_k=x_{k}
% p(x_{k}-Ax_{k-1})=exp(-.5(x_{k}-Ax_{k-1})^2./(Ddeltat+crlbest^2))
% pxx = 1-Pd; % also multiply with gaussian (p(x_{k}-Ax_{k-1})) but we implement that as a convolution.   
    %P(T<=delta t) = 1- exp(-koff deltat)
    if nargin < 4 || isempty(time)
        time = ones(size(LA1,1),1);
    end
       
    Drep = repmat(permute(D,[1 2]),[size(LA1,1) 1]); 
    timeRep = repmat(time,[1  probDim]);   
    
        
    pv = sum(dist./(2.*(2.*Drep.*timeRep+1.*(LA1.^2+LA2.^2))+eps)+0.5.*log(2.*pi.*(2.*Drep.*timeRep+1.*(LA1.^2+LA2.^2))+eps),2);
    if ~classic
        PdPv = pv-log(1-exp(-1.*kon).*kon);
    else
        PdPv = pv-log(1-koff);
    end
    PdPv(PdPv == 0) = 1e-10;
    
    PdPv(sqrt(sum(dist,2)) >3 | time > 3) = 1e6;
    %%
    PdPv(isinf(PdPv)) = 0;
    if isscalar(cc)
        pTrans = reshape(PdPv,[rr,cc]);
    else
        pTrans = sparse(rr,cc,PdPv);
    end
end



% function [ pTrans ] = transP(row,col,kon,koff,D,time,density,dist,probDim,LA1,LA2,rr,cc)
% 
% % function P supplies the transition probabilities as:
% % P(1,1)  s_{k-1}=0 s_k=0
% % P(1,2)  s_{k-1}=0 s_k=x_{k}
% % P(2,1)  s_{k-1}=x_{k-1} s_k=0
% % P(2,2)  s_{k-1}=x_{k-1} s_k=x_{k}
% % Created by Carlas Smith June 2014 (UMASS/TU-DELFT)
% 
% % TODO: add two sided boudndary crossing probaility for out of focus probability!!
% % Karlin, S. and Taylor H. M. (1975). A First Course in Stochastic Processes, 2nd ed., Academic
% % Press.
% % Jacod, J. and Protter P. (2004). Probability Essentials. 2nd edition. Springer. ISBN 3-540-
% % 43871-8.
% % P(|W_t| > b) = 4 \infsum_k (-1)^{k} phi((2k-1)b/sqrt(t)) assuming we at
% % zero.
% % function p = pwt(b,t)
% %     phi=@(x).5.*(1 + erf(x./sqrt(2)));
% %     k = 1:100*2*D*t;
% %     p = 4.*sum((-1).^(k).*phi((2.*k-1).*b./sqrt(2*D*t)));
% % end
% 
% if nargin < 8
%     dist = 1;
% end
% if nargin < 9
%     probDim = 2;
% end
% if nargin < 7 || nargin < 7
%     LA1 = 0; LA2 = 0; 
% end
% 
% if nargin < 12 || isempty(rr) ||  isempty(cc)
%     cc = 1; rr = 1;
% end
% % true = transition probabilities calculated by Pat Cutler.
% classic = false;
% 
% if row == 1 & col == 1
% if nargin < 4 || isempty(time)
%     time = ones(size(LA2,1),1);
% end
% % s_{k-1}=0 s_k=0
%     lr = sum(dist,2);
%     lr = (-log(1-density.*exp((time).*koff).*koff.*kon*exp(kon))).*(lr>0); % minimize costs of lr matrix
%     pTrans = lr;   
%     if isscalar(cc)
%         pTrans = reshape(pTrans,[rr,cc])';
%     else
%         pTrans = sparse(rr,cc,pTrans)';
%     end
% elseif row == 1 & col == 2
% % s_{k-1}=0 s_k=x_k
%         if nargin < 4 || isempty(time)
%             time = ones(size(LA2,1),1);
%         end
%         if ~classic
%             Pb = (-log(density)+(time).*koff-log(koff)+kon-log(kon));
%         else
%             Pb = (-log(density)-log(kon)-(time).*log(1-kon));
%         end
%         birthCost = Pb;
%         birthBlock = birthCost; %lower left
%         birthBlock(birthBlock==0) = NaN;
%         pTrans = birthBlock;
%         pTrans = diag(pTrans);
%         if ~isscalar(cc)
%             pTrans = sparse(pTrans);
%         end        
% elseif row == 2 & col == 1
% % s_{k-1}=x_{k-1} s_k=0
%         if nargin < 4 || isempty(time)
%             time = ones(size(LA1,1),1);
%         end
%         
%         %P(T<=delta t) = 1- exp(-koff deltat)
%         if ~classic
%             Pd = eps-real(log(density)+(time).*kon-log(kon)+koff-log(koff));
%         else
%             Pd = eps-(time).*log(1-kon);
%         end
%         deathCost = Pd ;
%         %generate upper right and lower left block
%         deathBlock = deathCost; %upper right
%         deathBlock(deathBlock==0) = NaN;
%         pTrans = deathBlock;
%         pTrans = diag(pTrans);
%         if ~isscalar(cc)
%             pTrans = sparse(pTrans);
%         end
% elseif row == 2 & col == 2
% % s_{k-1}=x_{k-1} s_k=x_{k}
% % p(x_{k}-Ax_{k-1})=exp(-.5(x_{k}-Ax_{k-1})^2./(Ddeltat+crlbest^2))
% % pxx = 1-Pd; % also multiply with gaussian (p(x_{k}-Ax_{k-1})) but we implement that as a convolution.   
%     %P(T<=delta t) = 1- exp(-koff deltat)
%     if nargin < 4 || isempty(time)
%         time = ones(size(LA1,1),1);
%     end
%        
%     Drep = repmat(permute(D,[1 2]),[size(LA1,1) 1]); 
%     timeRep = repmat(time,[1  probDim]);   
%     
%         
%     pv = sum(dist./(2.*(2.*Drep.*timeRep+(LA1.^2+LA2.^2))+eps)+0.5.*log(2.*pi.*(2.*Drep.*timeRep+(LA1.^2+LA2.^2))+eps),2);
%     if ~classic
%         PdPv = pv -real(log(1-(eps-real(log(density)+(time).*kon-log(kon)+koff-log(koff)))));   
%     else
%         PdPv = pv-log(1-koff);
%     end
%     PdPv(PdPv == 0) = 1e-10;
%     
%     PdPv(sum(dist,2) > 8) = 1e6;
%     %%
%     PdPv(isinf(PdPv)) = 0;
%     if isscalar(cc)
%         pTrans = reshape(PdPv,[rr,cc]);
%     else
%         pTrans = sparse(rr,cc,PdPv);
%     end
% end
% 
% end