function [ h ] = plotTracksC( trackLogA,closedNucleus, markers,lowresolution)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    if nargin < 4
        lowResolution = false;
    end
    numOfTracks=0;
    for i=1:length(trackLogA)
        numOfTracks=numOfTracks+length(trackLogA{i}.data);
    end
    
    co = jet(numOfTracks);

    h = figure;
    for i=1:length(trackLogA)
        plot3(zeros(3,1),zeros(3,1),zeros(3,1),'color','black','Marker', markers(i),'linewidth',2,'MarkerSize',3.75)
        hold on
        legendText{i} = ['(' num2str(i) ') '  trackLogA{i}.name];
    end
    atTrack = 1;
    for m=1:length(trackLogA)
         for i=1:length( trackLogA{m}.data)
                A = trackLogA{m}.data{i};
                plot3(A(:,3),A(:,2),A(:,1),'color',co(atTrack,:),'Marker', markers(m), 'MarkerSize',3.75)
                hold on
                atTrack=    atTrack+1;
         end
    end
    
    for i=1:size(closedNucleus,3)
%     [r,c,v] = ind2sub(size(squeeze(closedNucleusEdge(:,:,i-1))),find(permute(squeeze(closedNucleusEdge(:,:,i-1)),[1 2]) ==1 ));
%     plot(r,c,'or','MarkerSize',3.75)
%     hold on
        p(i) = patch(isosurface(repmat(logical(permute(squeeze(closedNucleus(:,:,i-1)),[2 1])),[1 1 500]),0.5));
        if lowresolution
            reducepatch(p(i),0.1)
        end
        set(p(i),'FaceColor',co(i,:),'EdgeColor','none','FaceAlpha',0.1);
        hold on
    end

    grid on
    legend(legendText)
    axis('tight')
end

