function [ trackLog trackLogA trackLogRaw diptrack] = importTrackFile( FileNameA,PathName, pixelSize )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    numOfTracks=length(FileNameA);

   for m=1:numOfTracks
        aa = importdata([PathName{m} FileNameA{m}]);
        a = str2double(aa.textdata(2:end,[3 4 9 2 1]));
        a(isnan(a)) = -1;
        minTrack = a(1,5);
        maxTrack = a(end,5);
        diptrack = zeros([maxTrack,max(a(:,4)),2]);
        
        trackLogRaw{m}.x = a(:,1);
        trackLogRaw{m}.y = a(:,2);
        trackLogRaw{m}.I = a(:,3);
        trackLogRaw{m}.t = a(:,4);
        trackLog{m}.data=[];
        for i=round(minTrack):round(maxTrack)
            trackLog{m}.name=FileNameA{m};
            
            trackLogRaw{i,m}.name=FileNameA{m};
            trackLogA{i,m}.x = [];
            trackLogA{i,m}.y = [];
            trackLogA{i,m}.I = [];
            trackLogA{i,m}.t = [];
            
            oneTrack = [];
            for j=1:size(a,1)
                trackNum = a(j,5);
                  if trackNum == i
                     diptrack(i,a(j,4)+1,:)= [a(j,1) a(j,2)]./pixelSize;
                     trackLogA{i,m}.x(end+1) = a(j,1)./pixelSize;
                     trackLogA{i,m}.y(end+1) = a(j,2)./pixelSize;
                     trackLogA{i,m}.I(end+1) = a(j,3);
                     trackLogA{i,m}.t(end+1) = a(j,4)+1;
                     dataone = aa.textdata(j+1,2:4);
                     oneTrack(end+1,1) = str2double(dataone{1});
                     oneTrack(end,2) = str2double(dataone{2})./pixelSize;
                     oneTrack(end,3) = str2double(dataone{3})./pixelSize;
                  end
            end
            trackLog{m}.data{end+1} = oneTrack;
        end
    end
end

