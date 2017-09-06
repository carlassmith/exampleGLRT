clear all
close all

load('ROIStack.mat')
PSFSigma=2;
Iterations=10;
[P, C, LL] = gpuGaussMLEv3(ROIStack,single(PSFSigma),Iterations,1,true ); 
P=P';
C=C';

[P1	, C1, LL] = gpuGaussMLEv3(ROIStack,single(PSFSigma),Iterations,4,true); 

%replace sigma info with fit results from fitType 2
P1=P1';    
C1=C1';
% P(:,6:7) = P1(:,6:6);
%%
figure;
plot(P1(:,5))
% C(:,5:6) = C1(:,5:6);

