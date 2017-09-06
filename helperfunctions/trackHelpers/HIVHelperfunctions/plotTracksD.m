function [ h ] = plotTracksD( trackLogA, markers,co,  IdxNucleus)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    
    if nargin < 4   
        IdxNucleus =-1;
    end
    
    numOfTracks=0;
    for i=1:length(trackLogA)
        numOfTracks=numOfTracks+length(trackLogA{i}.data);
    end
    
    h = figure;
    if IdxNucleus == -1
        subplot(1,2,1)
    end
    for i=1:size(trackLogA,2)
            plot(zeros(3,1),zeros(3,1),'color','black','Marker', markers(i),'LineWidth',2,'MarkerSize',3.75)
            hold on
            legendText{i} = ['(' num2str(i) ') '  trackLogA{i}.name];
    end
    grid on;
    
    atTrack = 1;
    alldistance = [];
    for m=1:size(trackLogA,2)
         for i=1:length( trackLogA{m}.distanceSM)
            if ~isempty(trackLogA{m}.data{i})                        
                 if IdxNucleus == -1
                     [val Idxmin] = min(min(trackLogA{m}.distanceSM{i}));
                     distance = trackLogA{m}.distanceSM{i}(:,Idxmin);
                 else
                    distance = trackLogA{m}.distanceSM{i}(:,IdxNucleus);
                 end
                    timeF = trackLogA{m}.data{i}(:,1);
                    plot(timeF,distance,'color',co(atTrack,:),'Marker', markers(m),'LineWidth',2,'MarkerSize',3.75)
                    hold on
                    atTrack=atTrack+1;
                alldistance=[alldistance;distance];
            end
        end
    end
    if IdxNucleus == -1 
        subplot(1,2,2)
        hist(alldistance,max(round(min(alldistance)+max(alldistance)),25));
        grid on;
    end
    legend(legendText)
    axis('tight')
end
