function [simcount, validcount] = CountEQ_In_DegBins(Allcat,validcat,pixlev,Mt)
%% remove all events smaller than Mt from the simulations and from the validcat
tol = 0.001;
Allcat(Allcat(:,4) < (Mt-tol),:) = [];

%%
binid = floor([Allcat(:,2)/pixlev(1), Allcat(:,3)/pixlev(2)]) + 1;
%%

Allcat(:,6:7)  = binid;
Allcat       = sortrows(Allcat,[6,7,5]);
[~,iafSim]   = unique(Allcat(:,6:7),'rows');
[usimtrid,ialSim] = unique(Allcat(:,6:7),'rows','last');


Prob2Occur        = cell(size(usimtrid,1),1);
for i=1:size(usimtrid,1)
    %i
    
    selsimid      = Allcat(iafSim(i):ialSim(i),5);
    [usimid,iaf]  = unique(selsimid);
    [usimid,ial]  = unique(selsimid,'last');
    count         = sort(ial-iaf+1);
    [countid,iaf] = unique(count);
    [countid,ial] = unique(count,'last');
    countid(:,2)  = ial-iaf+1;
    countid(:,3:4)  = repmat(usimtrid(i,:),[length(countid(:,1)),1]);
    Prob2Occur{i} = countid;
    check         = 1;
    
end

simcount = cell2mat(Prob2Occur);
%% assign each of the validation earthquake to triangle
if ~isempty(validcat)
    binid = floor([validcat(:,2)/pixlev(1), validcat(:,3)/pixlev(2)]) + 1;
    validcat(:,5:6)       = binid;
    validcat            = sortrows(validcat,[5,6]);
    [~,iafreal]         = unique(validcat(:,5:6),'rows');
    [urealtrid,ialreal] = unique(validcat(:,5:6),'rows','last');
    count               = ialreal-iafreal+1;
    validcount          = [urealtrid,count];
else
    validcount = [];
end
