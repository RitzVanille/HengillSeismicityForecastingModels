function [Allcat] = SimETAS_VContri(model,MomCat,Twin,safelim)
IP    = model.IP;
MomCat = [MomCat,IP];
if model.KDE_BKG == 1
    Adaptive_BW = repmat(10^model.partable.SMpars(1),size(model.IP));
else
    Adaptive_BW  =  model.smrad.^2;
    
end
b = model.partable.betaGlob/log(10);
enumbkg = sum(IP) * range(Twin) / range(MomCat(:,1));
nevt = poissrnd(enumbkg);

if nevt~=0
    % create the pool of the background earthquakes
    rn = rand(length(MomCat(:,1)),1);
    MomCat(IP < rn,:) = [];
    Adaptive_BW(IP < rn)    = [];
    % select from the random pool
    selind = randi(size(MomCat,1),nevt,1);
    selloc = MomCat(selind,2:3);
    selrad2 = Adaptive_BW(selind);
    % time and index of the background events
    nMomCat      = zeros(nevt,5);
    nMomCat(:,1) = rand(nevt,1)*(Twin(2)-Twin(1))+Twin(1);
    nMomCat      = sortrows(nMomCat,1);
    nMomCat(:,5) = [1:length(nMomCat(:,1))]';
    
    %% location of the background events
    randnum_loc   = rand(size(nMomCat,1),1);
    randnum_theta = rand(size(nMomCat,1),1)*360;
    rho = model.partable.SMpars(2);
    r             = sqrt(selrad2.*(1./(1-randnum_loc).^(1/rho)-1));
    degdist       = r*180/(6371*pi);
    [latout,lonout] = reckon(selloc(:,2),selloc(:,1),...
        degdist,randnum_theta);
    nMomCat(:,2)  =  lonout;
    nMomCat(:,3)  =  latout;
    %% simulate magnitude of the background events
    
    nbval = repmat(b,[nevt,1]);
    [nMomCat(:,4)] = GR_truncated_simulator_ver2(1,repmat(model.Mmax,[nevt,1]),...
        repmat(model.partable.Mpri - model.partable.dm/2,[nevt,1]),nbval);
    
    [nMomCat(:,4)] = BinMags(nMomCat(:,4),0,model.partable.dm);
    %%
    MomCat       =     nMomCat(:,1:4);
    Allcat       =     MomCat;
    Allcat(:,6)  =     0;
    gen          =     1;
    while size(MomCat,1) > 0
        
        [DaughtCat] = DirAft_ETAS(MomCat,model,Twin);
        if ~isempty(DaughtCat)
            DaughtCat(:,6) = gen;
            Allcat  =  [Allcat;DaughtCat];
            
        end
        MomCat  =  DaughtCat;
        gen     =  gen+1;
        
        if size(Allcat,1) > safelim
            disp('Safelimit exceeded; increase it or check the branching ratio')
            flg = 1;
            break
        end
    end
    
else
    Allcat  =  [];
end
end