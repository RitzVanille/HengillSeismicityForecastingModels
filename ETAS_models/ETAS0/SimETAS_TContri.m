function [Allcat] = SimETAS_TContri(model,MomCat,Twin,safelim)

%%
Allcat = [];
%%
gen = 1;
while size(MomCat,1) > 0
    [DaughtCat] = DirAft_ETAS(MomCat,model,Twin);
    if ~isempty(DaughtCat)
        DaughtCat(:,6) = gen;
        Allcat = [Allcat;DaughtCat];
    end
    size(DaughtCat);
    MomCat = DaughtCat;
    gen    = gen+1;
    if size(Allcat,1) > safelim
        disp('Safelimit exceeded; increase it or check the branching ratio')
        flg = 1;
        break
    end
    
end


end