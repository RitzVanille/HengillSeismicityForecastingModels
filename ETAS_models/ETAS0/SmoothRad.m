function [smrad]=SmoothRad(srcloc,nn,ip,minrad,jobind)

smrad=zeros(length(srcloc(:,1)),1);
temp=zeros(size(ip,1),3);
if isempty(jobind)
    jobind=[1,length(ip)];
end
for i=jobind(1):jobind(2)
%     if mod(i,100)==0
%         i;
%     end
    
    dist=(HaverSine_fast_ver2(srcloc(:,1:2),srcloc(i,1:2)))'; % distance in km
    temp(:,1:2)=[dist,ip];
    temp=sortrows(temp,1);
    temp(:,3)=cumsum(temp(:,2))-temp(1,2);
    %temp1=repmat((temp(:,3)),size(nn));
    [val,ind]=min(abs(temp(:,3)-nn));
    smrad(i)=temp(max(ind,2),1);
    
    
    check=1;
end
smrad(smrad<minrad)=minrad;

end