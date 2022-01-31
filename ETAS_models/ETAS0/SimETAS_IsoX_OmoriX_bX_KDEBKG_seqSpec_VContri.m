function [Allcat] = SimETAS_IsoX_OmoriX_bX_KDEBKG_seqSpec_VContri(model,MomCat,nsim,Twin,safelim)
tol = 0.001;
model.nprod = 0;
model.simtyp = 1;
%% select the background earthquaakes for simulation
IP    = model.IP;
smrad = model.smrad;
MomCat = [MomCat,IP];
MomCatold = MomCat;
if model.KDE_BKG == 1
    if model.SMkernel == 2
        Adaptive_BW = model.Adaptive_BW;
        mind = model.SMpars(3);
        Adaptive_BW(Adaptive_BW < mind^2) = mind^2;
    elseif model.SMkernel == 1
        Adaptive_BW = repmat(10^model.SMpars(1),size(model.IP));
    end
    check = 1;
elseif model.KDE_BKG == 2
    Adaptive_BW = model.Adaptive_BW;
    chekc = 1;
else
    
    Adaptive_BW  =  smrad.^2;
end

%% simulate magnitude
switch model.btyp
    case 1
        b = model.partable.betaGlob(1)/log(10);
        
    case 2
        b = model.partable.beta_bkg(1)/log(10);
end
enumbkg = sum(IP) * range(Twin) / range(MomCat(:,1));

%%
flg = 0;
i   = 1;
oldsmrad = smrad;
oldAdaptive_BW = Adaptive_BW;
while i <= nsim
    
    nevt = poissrnd(enumbkg);
    
    if nevt~=0
        % create the pool of the background earthquakes
        rn = rand(length(MomCat(:,1)),1);
        smrad(IP < rn)    = [];
        MomCat(IP < rn,:) = [];
        Adaptive_BW(IP < rn)    = [];
        % select from the random pool
        selind = randi(size(MomCat,1),nevt,1);
        selloc = MomCat(selind,2:3);
        selrad = smrad(selind);
        selrad2 = Adaptive_BW(selind);
        % time and index of the background events
        nMomCat      = zeros(nevt,5);
        nMomCat(:,1) = rand(nevt,1)*(Twin(2)-Twin(1))+Twin(1);
        nMomCat      = sortrows(nMomCat,1);
        nMomCat(:,5) = [1:length(nMomCat(:,1))]';
        
        if model.KDE_BKG == 0
            % generate a distance and azimuth
            sigma  =  repmat([1,0;0,1],[1,1,nevt]);
            for i  =  1:nevt
                sigma(:,:,i)  =  sigma(:,:,i)*(selrad(i)^2);
            end
            [tempx]  =  mvnrnd(zeros(nevt,2),sigma);
            R        =  sqrt(sum(tempx.^2,2));
            degdist  =  R*180/(6371*pi);
            theta    =  rand(size(R))*360;
            [latout,lonout]  =  reckon(selloc(:,2),selloc(:,1),degdist,theta);
            
            nMomCat(:,2)  =  lonout;
            nMomCat(:,3)  =  latout;
        else
            randnum_loc   = rand(size(nMomCat,1),1);
            randnum_theta = rand(size(nMomCat,1),1)*360;
            rho = model.SMpars(2);
            r             = sqrt(selrad2.*(1./(1-randnum_loc).^(1/rho)-1));
            degdist       = r*180/(6371*pi);
            [latout,lonout] = reckon(selloc(:,2),selloc(:,1),...
                degdist,randnum_theta);
            nMomCat(:,2)  =  lonout;
            nMomCat(:,3)  =  latout;
            check = 1;
        end
        
        
        %% simulate magnitude
        if model.btyp == 3.1
            check = 1;
            [globvars] = GlobvarsDistSSM_ver4(model.bmodel.pivots,[lonout,latout],model.bmodel.maxdist);
            Cglob = [1:size(lonout,1)]';
            iaf = (Cglob-1)*model.bmodel.maxdist + 1;
            ial = Cglob*model.bmodel.maxdist;
            nbval = zeros(size(Cglob));
            
            betas = [model.bmodel.pars(2:end-1)';repmat(model.bmodel.pars(end),[model.bmodel.Npivots2,1])];
            lambda = model.bmodel.pars(1);
            for i = 1:length(Cglob)
                selbeta = betas(globvars(iaf(i):ial(i),3));
                wts = exp(-lambda*globvars(iaf(i):ial(i),1));
                wts = wts/sum(wts);
                nbval(i) = sum(wts.*selbeta);
            end
            
            nbval = nbval/log(10);
            [nMomCat(:,4)] = GR_truncated_simulator_ver2(1,repmat(model.Mmax(1),[nevt,1]),...
                repmat(model.M0-model.partable.dm(1)/2,[nevt,1]),nbval);
            
        else
            nbval = repmat(b,[nevt,1]);
            [nMomCat(:,4)] = GR_truncated_simulator_ver2(1,repmat(model.Mmax(1),[nevt,1]),...
                repmat(model.M0 - model.partable.dm(1)/2,[nevt,1]),nbval);
            
            
        end
        [nMomCat(:,4)] = BinMags(nMomCat(:,4),0,model.dm);
        %%
        MomCat       =     nMomCat(:,1:4);
        MomCat(:,6)  =     0;
        MomCat(:,7)  =     1;
        MomCat(:,8)  = [1:size(MomCat,1)]';
        Allcat       =     MomCat;
        gen          =     1;
        model.EQfam = MomCat(:,7);
        while size(MomCat,1) > 0
            %%
            [DaughtCat] = DirAft_ETAS_IsoX_OmoriX_bX_KDEBKG_seqSpec(MomCat,model,Twin);
            if ~isempty(DaughtCat)
                DaughtCat(:,6) = gen;
                Allcat  =  [Allcat;DaughtCat];
                if model.seqSpec == 1
                    model.EQfam = DaughtCat(:,7);
                end
            end
            MomCat  =  DaughtCat;
            gen     =  gen+1;
            
            if size(Allcat,1)  >  safelim
                flg  =  1;
                break
            end
        end
        
        if flg  ==  0
            if ~isempty(Allcat)
                Allcat(Allcat(:,4) < (model.Mt-tol),:)  =  [];
                %Allcat                            =   sortrows(Allcat,1);
            end
            i  =  i+1;
            
        end
        MomCat  =  MomCatold;
        smrad   =  oldsmrad;
        Adaptive_BW   =  oldAdaptive_BW;
        flg     =  0;
    else
        Allcat  =  [];
    end
    
    
end

end