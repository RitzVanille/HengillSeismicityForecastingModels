function [globvars,Globcell]=GlobVars(catsrc,catrec,srclen,discLatLon)
%% this function computes the distance between pair of events in the source and
% the target catalog

% The maximum allowed distance is defined in terms of number of source
% lengths

% one source length of the event is defined using the wells and coppesmith
% equations

% distances between the causal pairs are stored

% the indexes of the events from the source and target catalog are stored

%%
%{
globvars:
1. distanceij^2
2. tij
3. magnitude of the source
4. j: index of the target
5. i: index of the source
6. distance in number of source length
7. magnitude of the target


%}

%%
Globcell=cell(length(catsrc(:,1)),1);
for i = 1:length(catsrc(:,1))
    
    ind=catrec(:,1)>catsrc(i,1);
    if sum(ind)>0
        dist=(HaverSine_fast_ver2(catrec(ind,2:3),catsrc(i,2:3)))';
        tmdif=catrec(ind,1)-catsrc(i,1);
        magrec=catrec(ind,4);
        [~,SSRL] = coppersmith_ver2(catsrc(i,4),4);
        
        Globcell{i,1}=[dist.^2,tmdif,catsrc(i,4)*ones(length(dist),1),...
            find(ind),i*ones(length(dist),1),(dist)/SSRL,magrec];
        Globcell{i,1}(Globcell{i,1}(:,6)>srclen,:)=[];
    else
        Globcell{i,1}=[];
        
    end
    %toc
    check=1;
end
globvars=cell2mat(Globcell);
globvars=sortrows(globvars,4);

if discLatLon == 1
    globvars(globvars(:,1) == 0,:) = []; % when the catalog is imprecise, discrete lat long leads to many 0 distances.
end

end

