function [globvars,maxrad] = GlobvarsDistSSM_ver3(loc,nsample)
%% this function computes the distance between pair of events 
% this input is required for the inversion of parameters of the smoothed
% seismicity model

% % To avoid very large matrix, only n nearest neighbours are stored for each
% event

%{
inputs
----------------------
loc: 
1. longitude
2. latitude

nsamle:
number of nearest neighbours, choose what you can afford

outputs
----------------------

globvars:
1. distanceij^2
2. index i
3. index j

maxrad:
maximum influence radius of each event

%}

%%
allids = [1:length(loc(:,1))]';
Globcell=cell(length(loc(:,1)),1);
dispiter = round(length(loc(:,1))/100);

for i = 1:size(loc,1)
   
    dist = (HaverSine_fast_ver2(loc,loc(i,:)))';
    
    Globcell{i,1} = sortrows([dist.^2,i*ones(length(dist),1),allids],1);
    Globcell{i,1} = Globcell{i,1}(1:nsample,:);
    check = 1;
end

globvars=cell2mat(Globcell);
globvars=sortrows(globvars,3);
globvars(globvars(:,1) == 0, : ) = [];


globvars1 = sortrows(globvars,2);
[Cglob1,iaf1]  = unique(globvars1(:,2));
[Cglob1,ial1]  = unique(globvars1(:,2),'last');
maxrad = zeros(size(Cglob1));
for i1 = 1:length(Cglob1)
    maxrad(i1) = max(globvars1(iaf1(i1):ial1(i1),1));
end
check = 1;
end