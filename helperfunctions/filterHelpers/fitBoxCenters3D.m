function [ RawFitResults ] = fitBoxCenters3D(data,boxSize,boxCenters,sigmaPSF,maxCudaFits,iterations,fitType)

%% Extract subregions for fits

% params{1}.sigma_xy=1.9965;
% params{1}.sigma_z=2.9588;
% 
% params{1}.Npixels_xy = 3*(2*params{1}.sigma_xy+1);
% params{1}.Npixels_z = 3*(2*params{1}.sigma_z+1);

% [ROIStack XRoiStart YRoiStart] = cMakeSubregions(boxCenters(:,1),boxCenters(:,2),boxCenters(:,3),boxSize,single(permute(data,[2 1 3])));
[ROIStack, top,left,depth] = cMakeSubregionsv2(boxCenters(:,1),boxCenters(:,2),boxCenters(:,3),ceil(boxSize(1)),single(data),ceil(boxSize(2)));
%[ROIStack, top,left,depth] = cMakeSubregionsv2(boxCenters(:,2),boxCenters(:,1),boxCenters(:,3),ceil(boxSize(1)),single(data),ceil(boxSize(2)));

XRoiStart = single(top);
YRoiStart = single(left);
ZRoiStart = single(depth);

NROI=size(ROIStack,4); %number of ROIs

%fit data (fitType ?)
[P, C, L] = gpuGaussMLEv3(permute(single(ROIStack),[2 1 3 4]),single(sigmaPSF),iterations,3,true,[maxCudaFits]); 
% [P, C, L]=gaussmlev2(ROIStack(:,:,st:en),sigmaPSF,iterations,1);
%fit data (fitType ? estimate sigma)
[P1, C1, ~] = gpuGaussMLEv3(permute(single(ROIStack),[2 1 3 4]),single(sigmaPSF/2),iterations,5,true,[maxCudaFits]); 
% [P1, C1]=gaussmlev2(ROIStack(:,:,st:en),sigmaPSF,iterations,2);
%replace sigma info with fit results from fitType 2
P(6,:) = P1(6,:);
P(7,:) = P1(7,:);
C(6,:) = C1(6,:);
C(7,:) = C1(7,:);

P=P';
C=C';
L=L';
%Populate RawFitResults
significance_default = 0.05;

%ROI starting positions
RawFitResults.RoiStart = [XRoiStart YRoiStart  ZRoiStart];
RawFitResults.ROIStack = ROIStack;

%estimted parameters
RawFitResults.Coord(:,1)=double(P(:,1))+double(XRoiStart);
RawFitResults.Coord(:,2)=double(P(:,2))+double(YRoiStart);
RawFitResults.Coord(:,3)=double(P(:,3))+double(ZRoiStart);
RawFitResults.Photons=double(P(:,4));
RawFitResults.Bg=double(P(:,5));
RawFitResults.SigmaXY=double(P(:,6));
RawFitResults.SigmaZ=double(P(:,7));

%CRB for parameters
RawFitResults.CRLB_STD(:,1)=double(C(:,1));
RawFitResults.CRLB_STD(:,2)=double(C(:,2));
RawFitResults.CRLB_STD(:,3)=double(C(:,3));
RawFitResults.Photons_STD=double(C(:,4));
RawFitResults.Bg_STD=double(C(:,5));
RawFitResults.SigmaXY_STD=double(C(:,6));
RawFitResults.SigmaZ_STD=double(C(:,7));

[ ~,~,pfa_adj ]=fdr_bh(reshape(double(L(:,3)),[prod(size(L(:,3))) 1]),significance_default,'dep','no');

RawFitResults.PFA = min(pfa_adj,1);
RawFitResults.LL=double(L);
