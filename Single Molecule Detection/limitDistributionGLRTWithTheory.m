% Code for the publication:
%  C.S. Smith, S. Stallinga, K.A. Lidke, B. Rieger, D. Grünwald, 
%  Probability-based particle detection that enables threshold-free and robust in vivo single molecule tracking,
%  Molecular Biology of the Cell, 2015 accepted. mbc.E15-06-0448
%
% Copyright 2015, University of Massachusetts, USA
%
% Carlas Smith, October 2015
%
% This code requires the toolbox DIPimage and CUDA 6.5 to be installed
% The windows watchdog that will stop the GPU execution if not disabled! 
% (TdrLevel=0,TdrDelay=180; https://msdn.microsoft.com/en-us/library/windows/hardware/ff569918(v=vs.85).aspx) 
% (http://www.diplib.org/download; https://developer.nvidia.com/cuda-toolkit-65)

close all
clearvars

addpath(genpath('helperfunctions'))
addpath(genpath('../helperfunctions/'))


%% simulation parameters
bgc=[10 15 20]; %background fluorescence in photons/pixel/frame
I0 =[150 100 75 50]; %expected photons/frame

sigma=1.39; %PSF sigma in pixels
Npixels=ceil(3*(2*sigma+1)); %linear size of fit region in pixels. 

% Make sure we have enough pixels
part = Npixels/2-floor(Npixels/2);

if part == 0;
    Npixels=Npixels+1;
end

Nsamples = 1024;
Nsims = 16;
resolution = 500;
%% calculate non-centrality parameters for non central chi square distribution
for ii=1:length(I0);
    for bb=1:length(bgc)
        [crlb I] = fisherInformationGauss(Npixels,sigma,I0(ii),bgc(bb));
         lambda(ii,bb) = (I0(ii)).^2.*1./crlb(3,3); 
    end
end

%% Calculate theoretical limit distribution and do ratio test on simulated data
PFA = logspace(-3,0,20);

for ii=1:length(I0);
      % Calculate theoretcical probability of detection for fixed false
      % positive probabilty
      for bb=1:length(bgc)
        QinvPFA(ii,bb,:) = norminv(1-PFA/2,0,1);
        PD(ii,bb,:) = (1-normcdf(QinvPFA(ii,bb,:),-sqrt(lambda(ii,bb)),1))+(1-normcdf(QinvPFA(ii,bb,:),sqrt(lambda(ii,bb)),1));
      end
      
      for q = 1:Nsims
          for bb=1:length(bgc)
              % Do ratio test on simulated data
              [pMapP, pMapFA] = doBgIRatio(Npixels,sigma,I0(ii),bgc(bb),Nsamples); 
              PFAp(ii,bb,:,q)=pMapFA;
              
              % Set fixed false postive rate
              for i=1:size(PFA,2)
                 PFAp(ii,bb,i)=PFA(i);
                 decisionclas(ii,bb,i,:,q) = pMapFA <= PFA(i);
            end
          end 
      end
end


%% Plot results and calculate ROC curve.
cm={'r','b','g','c','k'};
ln={'-','--','s-',':'};
h = figure;
plot(0,0,'k-','linewidth',1)
hold on
plot(0,0,'ko','linewidth',1)

for ii=1:length(I0);
    plot(0,0,[ln{ii} 'k'],'linewidth',1)
    hold on
end

decisionm = mean(decisionclas,5);
for ii=1:length(I0);
      for bb=1:length(bgc)
        hold on
        loglog(PFA,squeeze(PD(ii,bb,:)),[cm{bb} ln{ii}],'linewidth',2)
      end
      for bb=1:length(bgc)
            detectFalseStd(ii,bb,:)= sqrt(sum(decisionm(ii,bb,:,1:Nsamples).^2,4)./Nsamples);
            detectTrueStd(ii,bb,:)= sqrt(sum(decisionm(ii,bb,:,Nsamples+1:end).^2,4)./Nsamples);
            detectTrue(ii,bb,:) = sum(decisionm(ii,bb,:,Nsamples+1:end),4)./Nsamples;   
            detectFalse(ii,bb,:) = sum(decisionm(ii,bb,:,1:Nsamples),4)./Nsamples;   

            errorbar(squeeze(detectFalse(ii,bb,:)),squeeze(detectTrue(ii,bb,:)), squeeze(detectFalseStd(ii,bb,:)),[cm{bb} 'o'],'linewidth',2);
            hold on
      end 
end
ylabel('True positive rate','fontsize',16)
xlabel('False positive rate','fontsize',16)  

clear lengendtext;
lengendtext{1} = 'Theory';
lengendtext{2} = 'Simulation';

for ii=1:length(I0);
    lengendtext{ii+2} = sprintf('Intensity %d',round(I0(ii)));    
end
    
for bb=1:length(bgc)
    lengendtext{end+1} = sprintf('Background %d',round(bgc(bb)));    
end

handleL3 = legend(lengendtext,'Location','SouthEast');
set(handleL3,'FontSize',10);
grid on
set(gca,'xscale','log');
set(gca,'yscale','log');
xlim([10^-2 10^(-1)]);
ylim([10^-1 10^0]);

