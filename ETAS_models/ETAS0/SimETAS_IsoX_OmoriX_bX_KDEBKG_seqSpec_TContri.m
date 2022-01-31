function [Allcat] = SimETAS_IsoX_OmoriX_bX_KDEBKG_seqSpec_TContri(model,MomCat,nsim,Twin,safelim,nprod)
MomCatold = MomCat;
tol = 0.001;
i   = 1;
flg = 0;
if model.seqSpec == 1
    oldEqFam = model.EQfam;
end
while i <= nsim
    if model.seqSpec == 1
        
        model.EQfam = oldEqFam;
    end
    Allcat = [];
    %%
    gen = 1; SimTyp = 2;
    if model.isotyp == 2
        model.simtyp = 2;
    else
        model.simtyp = 1;
    end
    
    model.nprod = nprod;
    MomCat(:,8) = [1:size(MomCat,1)]';
    
    check = 1;
    while size(MomCat,1) > 0
        [DaughtCat] = DirAft_ETAS_IsoX_OmoriX_bX_KDEBKG_seqSpec(MomCat,model,Twin);
        
        if ~isempty(DaughtCat)
            DaughtCat(:,6) = gen;
            Allcat = [Allcat;DaughtCat];
            if model.seqSpec == 1
                model.EQfam = DaughtCat(:,7);
            end
        end
        size(DaughtCat);
        MomCat = DaughtCat;
        gen    = gen+1;
        if size(Allcat,1) > safelim
            flg = 1;
            break
        end
        model.simtyp = 1;
        model.nprod = 0;
        
        
        %toc
    end
    
    if flg == 0
        if ~isempty(Allcat)
            Allcat(Allcat(:,4) < (model.Mt-tol),:) = [];
            %Allcat   =  sortrows(Allcat,1);
        end
        i = i+1;
    end
    flg    = 0;
    MomCat = MomCatold;
end
check = 1;
end