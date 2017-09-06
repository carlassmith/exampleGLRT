function [ trackLogA trackLogRaw diptrack] = importTrackFile( FileNameA,PathName, pixelSize )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    numOfTracks=length(FileNameA);

   for m=1:numOfTracks
        a = importdata([PathName{m} FileNameA{m}]);
        minTrack = a(1,5);
        maxTrack = a(end,5);
        diptrack = zeros([maxTrack,max(a(:,4)),2]);
        
        trackLogRaw{m}.x = a(:,1);
        trackLogRaw{m}.y = a(:,2);
        trackLogRaw{m}.I = a(:,3);
        trackLogRaw{m}.t = a(:,4);

        for i=round(minTrack):round(maxTrack)
            trackLogA{i+1,m}.name=FileNameA{m};
            trackLogRaw{i+1,m}.name=FileNameA{m};
            trackLogA{i+1,m}.x = [];
            trackLogA{i+1,m}.y = [];
            trackLogA{i+1,m}.I = [];
            trackLogA{i+1,m}.t = [];
            
            oneTrack = [];
            for j=1:size(a,1)
                trackNum = a(j,5);
                  if trackNum == i
                     diptrack(i+1,a(j,4)+1,:)= [a(j,1) a(j,2)];
                     trackLogA{i+1,m}.x(end+1) = a(j,1);
                     trackLogA{i+1,m}.y(end+1) = a(j,2);
                     trackLogA{i+1,m}.I(end+1) = a(j,3);
                     trackLogA{i+1,m}.t(end+1) = a(j,4)+1;
                  end
            end
        end
    end
end

