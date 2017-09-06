function [ mpval, pTracks ] = pTransTracks( tracks, out, costMatF2F  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Good enough
D = costMatF2F.D;
kon = costMatF2F.kon;
koff = costMatF2F.koff;
density = costMatF2F.density;

for i=1:size(tracks,1)
    pTracks(i).pval=[];
    for j=1:size(out(i).Frame,2)-1
     
            deltat = out(i).Frame(j+1)-out(i).Frame(j);
            deltax = out(i).X(j+1)-out(i).X(j);
            deltay = out(i).Y(j+1)-out(i).Y(j);
            LA1 = [out(i).std_x(j) out(i).std_y(j)];
            LA2 = [out(i).std_x(j+1) out(i).std_y(j+1)];            
            pTracks(i).pval(j)= (transP(2,2,kon,koff,D,deltat,density,[deltax deltay],2,LA1,LA2)-log(out(i).pH1(j))-log(out(i).pH1(j+1)));
      
    end
    if size(out(i).Frame,2) == 1
        pTracks(i).pval(j)=-log(out(i).pH1(j));
    end
    if ~isempty(pTracks(i).pval)
%         mpval(i) = mean(exp(-pTracks(i).pval));
          mpval(i) = mean(out(i).pH1);
    else
        mpval(i)=0;  
    end
end

% The right way to do it:
% pTransTracks = zeros(size(tracks,1),1);
% for i=1:size(tracks,1)
%     deltaT = out(i).Frame;
%     dist = [diff(out(i).X) diff(out(i).Y)];
%     if out(i).Frame(1) > 1
% %           ~  -> on 1,2
%          disp(' -> on 1,2 ')
%         pTransTracks(i) = pTransTracks(i) + transP(1,2,kon,koff,D,1,size(boxCenters,1)./szP);
%     end
%     if size(out(i).Frame,2) > 1
%         for j=1:size(out(i).Frame,2)-1
%             if out(i).Frame(j) - out(i).Frame(j+1) > 1
%     %             off -> on 1,2 
%                 disp(' off -> on 1,2 ')
%                 pTransTracks(i) = pTransTracks(i) + transP(1,2,kon,koff,D,out(i).Frame(j) - out(i).Frame(j+1),size(boxCenters,1)./szP);
%             else
%     %               on-> on 2,2
%                 disp(' on -> on 2,2 ')
%                 LA1 = [out(i).std_x(j) out(i).std_y(j)];
%                 LA2 = [out(i).std_x(j+1) out(i).std_y(j+1)];
%                 pTransTracks(i) = pTransTracks(i) + transP(2,2,kon,koff,D,[],size(boxCenters,1)./szP,dist(j,:),2,LA1,LA2);
%             end
%         end   
%         if out(i).Frame(end) < frames
%     %         -> off 2,1
%                 disp(' off -> on 1,2 ')
%                 pTransTracks(i) = pTransTracks(i) + transP(2,1,kon,koff,D,out(i).Frame(end-1) - out(i).Frame(end),size(boxCenters,1)./szP);
%         end
%     end
%     val(i) = mean(out(i).pH1);
% end
pVal = -1;

end

