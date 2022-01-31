function [output] = Wrapper_Sim_ETAS_verCoSeismiq(nsim,maxtm,InvUpdateDur,SimUpdateDur,...
    Maux,Mpri,KDE_BKG,OmTyp,btyp,isotyp,prodtyp,srclen,mthres,nprod,seqSpec,minFamSz,tresol,jobindex)
%{

    papers in which these codes were developed

@article{nandan2019forecasting,
  title={Forecasting the rates of future aftershocks of all generations is essential to develop better earthquake forecast models},
  author={Nandan, Shyam and Ouillon, Guy and Sornette, Didier and Wiemer, Stefan},
  journal={Journal of Geophysical Research: Solid Earth},
  volume={124},
  number={8},
  pages={8404--8425},
  year={2019}
}

@article{nandan2021global,
  title={Global models for short-term earthquake forecasting and predictive skill assessment},
  author={Nandan, Shyam and Kamer, Yavor and Ouillon, Guy and Hiemer, Stefan and Sornette, Didier},
  journal={The European Physical Journal Special Topics},
  volume={230},
  number={1},
  pages={425--449},
  year={2021},
  publisher={Springer Berlin Heidelberg}
}

@article{nandan2021seismicity,
  title={Is seismicity operating at a critical point?},
  author={Nandan, Shyam and Ram, Sumit Kumar and Ouillon, Guy and Sornette, Didier},
  journal={Physical Review Letters},
  volume={126},
  number={12},
  pages={128501},
  year={2021},
  publisher={American Physical Society}
}
    
    
    
    [output] = Wrapper_Sim_ETAS_verCoSeismiq(100,inf,30,30,...
    0.3,0.3,1,1,1,1,0,50,0,0,0,0,30,12)
    
    %}
    %tmbegin = tic;
    %%
    
    
    
    baseoutfolder = './outputvars/';
    
    plt = 0;
    cursim = 1;
    penaltyFac = 0.4;
    nn         = 700;
    compresssyn = 0;
    
    usebr     = 0 ;
    applyBoundary = 1;
    eqArea    =   0;
    syncat    =   0;
    ssmresol = 0.1;
    %tresol    =   30;
    tol        =  10^-3;
    %savesim    =  1;
    onlysim   = 0;
    if syncat == 0
        safelim   =   200000;
    else
        safelim = 10^8;
    end
    
    
    catnam   =  'Coseismiq';
    Regionnam = 'Iceland-Hengill-Active_pval0.1';
    Prefix = strcat(catnam,'_',Regionnam);
    Mmax = 8.5;
    srclen0 = 50;
    begsttm = datenum(2020,2,1);
    begsttm_inv = datenum(2020,2,1);
    Sresol     =  [0.005];
    Mt         =  [0.3];
    savcat = 1;
    applyBoundary = 1;
    Mmaxlim = [7,8];
    
    cat0   =  importdata(strcat(Prefix,'.mat'));
    
    
    styraux = cat0.styraux;
    
    
    simfoldnam = strcat('./outputvars/',Prefix,'/Sim_updateDur',num2str(SimUpdateDur),'/');
    
    mkdir(simfoldnam)
    
    modpath  =   strcat(baseoutfolder,Prefix,'/Inv_updateDur',num2str(InvUpdateDur),'/');
    
    %% decide which parameter to use
    endsttm_inv = max(cat0.pricat.datenum);
    sttmlist_inv = [begsttm_inv:InvUpdateDur:endsttm_inv];
    sttmlist_sim = [begsttm_inv:SimUpdateDur:endsttm_inv];
    dur = tresol;
    
    njobindex = jobindex;
    
    
    sttmid = mod(jobindex,length(sttmlist_sim)) + 1;
    
    sttm_sim = sttmlist_sim(sttmid);
    
    [Idx,d] = knnsearch(sttmlist_inv',sttm_sim,'K',2);
    
    if sum(d==0) == 1
        sttmid_inv = Idx(1);
    else
        sttmid_inv = min(Idx);
    end
    sttm_inv = sttmlist_inv(sttmid_inv);
    sttm = sttm_sim;
    endtm = sttm + dur;
    
    
    rng(njobindex)
    %% deal with special cases
    invisotyp = isotyp;
    invbtyp = btyp;
    %% load the model
    modnam     =   strcat(modpath,'Pars_sttm',num2str(sttm_inv),'_Maux',num2str(Maux),...
        '_Mpri',num2str(Mpri),'_KDEBKG',num2str(KDE_BKG),'_OmTyp',num2str(OmTyp),'_btyp',num2str(invbtyp)...
        ,'_isotyp',num2str(invisotyp),'_srclen',num2str(srclen),'_mthres',num2str(mthres),...
        '_prodtyp',num2str(prodtyp),'_seqSpec',num2str(seqSpec),'_minFamSz',num2str(minFamSz),'.mat')
    
    
    
    model       =    importdata(strcat(modnam));
    model.isotyp  =   isotyp;
    model.btyp  =   btyp;
    
    
    disp(Prefix)
    disp(datestr(sttm_sim))
    %% prepare the training catalog
    dm = cat0.dm;
    pricat = [cat0.pricat.datenum,cat0.pricat.longitude,cat0.pricat.latitude,cat0.pricat.magnitude];
    [pricat(:,4)]  =  BinMags(pricat(:,4),0,dm);
    auxcat = [cat0.auxcat.datenum,cat0.auxcat.longitude,cat0.auxcat.latitude,cat0.auxcat.magnitude];
    [auxcat(:,4)]  =  BinMags(auxcat(:,4),0,dm);
    styrpri = cat0.styrpri;
    if styrpri<2100
        styrpri  =  datenum(styrpri,1,1);
    else
        styrpri = cat0.styrpri;
    end
    
    catsrc =  auxcat(auxcat(:,1) >= styraux & auxcat(:,1) < sttm & auxcat(:,4)>= Maux,:);
    catrec =  pricat(pricat(:,1) >= styrpri & pricat(:,1) < sttm_inv & pricat(:,4)>= Mpri,:);
    %catsrc(end+1,1:4) = [sttm-0.0001,8.5417,47.3769,5];
    %% validation catalog
    
    %for j = 1:length(tresol)
    tmpendtm = sttm + tresol;
    for i = 1:length(Mt)
        validcat{1,i}  =  auxcat(auxcat(:,4) >= Mt(i) & auxcat(:,1) >= sttm & auxcat(:,1) < tmpendtm,:); % select only events above auxilliary
    end
    %end
    %%
    
    IP    =   model.IP;
    IP(catrec(:,1) > sttm_inv,:) = [];
    model.IP = IP;
    smrad    = model.smrad;
    smrad(catrec(:,1) > sttm_inv,:) = [];
    model.smrad = smrad;
    Adaptive_BW    = model.Adaptive_BW;
    Adaptive_BW(catrec(:,1) > sttm_inv,:) = [];
    model.Adaptive_BW = Adaptive_BW;
    
    
    %%
    
    
    model.lambdai = [];
    
    
    model.Mmax = Mmax;
    model.Mmin = Maux;
    model.Mt   = min(Mt);
    model.M0   = Mpri;
    model.M1   = Maux;
    model.usebr = usebr;
    model.catsrc = catsrc;
    model.catrec = catrec;
    model.isotyp = isotyp;
    model.prodtyp = prodtyp;
    model.mthres = mthres;
    
    
    %% Model description
    switch seqSpec
        case 0
            SEQ_Text = 'Global Parameters' % space and temporally homogenous poisson process
            model.seqSpec = 0;
        case 1
            SEQ_Text = 'Sequence specific triggering parameters' % space varying temporally homogenous poisson process
    end
    
    %%
    switch KDE_BKG
        case 0
            KDE_Text = 'STHPP' % space and temporally homogenous poisson process
        case 1
            KDE_Text = 'SV_THPP' % space varying temporally homogenous poisson process
        case 2
            KDE_Text = 'SV_THPP with adaptive kernel' % space varying temporally homogenous poisson process
    end
    %%
    switch prodtyp
        case 0
            ProdText  =  'Parameteric' % temporary
        case 1
            ProdText = 'Bilinear Parametric with same spatial kernel'
            %isotyp = 1;
            mthres
        case 2
            ProdText = 'Bilinear Parametric with different spatial kernel'
            %isotyp = 1;
        case 3
            ProdText = 'Linear magnitude dependence of alpha'
        case 301
            ProdText = 'RELU magnitude dependence of alpha'
        case 302
            ProdText = 'TanH tapered RELU magnitude dependence of alpha'
        case 303
            ProdText = 'alpha(m) = a0 * exp(a1*m)'
        case 304
            ProdText = 'alpha(m) = a0 * exp(a1*m)'
    end
    
    %%
    switch OmTyp
        case 1
            OmText = 'Exponentially tapered Omori Model'
        case 2
            OmText = 'p(m) and c(m)'
        case 3
            OmText = 'p(m)'
        case 4
            OmText = 'p(m) c(m) and tau(m)'
        case 0
            OmText = 'Pure Omori'
            
    end
    %%
    switch btyp
        case 1
            bText = 'Same b for background and aftershocks' % temporary
        case 2
            bText = 'generalized verejones model' % temporary
    end
    
    %isotyp = 3;
    
    %%
    
    switch isotyp
        case 1
            isoText = 'Isotropic' % temporary
        case 2
            isoText = 'Anisotropic' % temporary
            
        case 3
            isoText = 'Isotropic 2.0' % temporary
            
        case 4
            isoText = 'Magnitude dependent rho' % temporary
    end
    
    %%
    disp('Starting Simulations...')
    lencat = zeros(nsim,length(Mt));
    % lencat_aft1 =  zeros(nsim,1);
    % lencat_aft2 =  zeros(nsim,1);
    
    Allcat = cell(nsim,1);
    tic
    for i = 1:nsim
        %tic
        %tic
        [aft1] = SimETAS_IsoX_OmoriX_bX_KDEBKG_seqSpec_TContri(model,catsrc,1,[sttm,endtm],safelim,nprod);
        %toc
        [aft2] = SimETAS_IsoX_OmoriX_bX_KDEBKG_seqSpec_VContri(model,catrec,1,[sttm,endtm],safelim);
        
        %     lencat_aft1(i,1) = size(aft1,1);
        %     lencat_aft2(i,1) = size(aft2,1);
        aft    = [aft1;aft2];
        if ~isempty(aft)
            %aft(:,1:4) = round(aft(:,1:4),2);
            aft(:,5)   =  i;
            Allcat{i}     =  aft(:,1:7);
        end
        tmelapsed      =  toc;
        nsimcomp       =  i;
        if mod(i,nsim/100) == 0
            %clc
            disp(strcat(['Simulation Stage:: ', ...
                num2str(round(i * 100/nsim)),'%']))
            toc
            %dispMmax = max(Allcat(:,4))
        end
        %for k = 1: length(tresol)
        %naft = aft(aft(:,1)<(sttm+tresol),:);
        %lencat{k} = zeros(nsim,length(Mt));
        for j = 1:length(Mt)
            lencat(i,j)  =  sum(aft(:,4) >= (Mt(j)-tol));
        end
        %end
        check      =  1;
        if tmelapsed>=maxtm
            break
        end
    end
    disp('Simulations over.')
    
    %Allcat{10001,1} = [];
    Allcat = cell2mat(Allcat);
    
    %% apply the in polygon function if Cregion is available
    if syncat == 1
        Allcat = sortrows(Allcat,1);
    end
    
    if applyBoundary == 1
        if ~isempty(cat0.Cregion)
            [IN] = inpolygon(Allcat(:,2),Allcat(:,3),cat0.Cregion(:,1),cat0.Cregion(:,2));
            Allcat(~IN,:) = [];
        end
        
    end
    %% count the number of earthquakes in each spatial bin
    Allcat(:,1) = Allcat(:,1) - sttm;
    
    for i = 1:length(Sresol)
        for j = 1:length(Mt)
            [SimCount{i,j}, ValidCount{i,j}] = CountEQ_In_DegBins(Allcat,validcat{1,j},[Sresol(i),Sresol(i)],Mt(j));
        end
    end
    %%
   
    output.pars       =      model.partable;
    if savcat == 1
        output.Allcat = Allcat;
    end
    output.sttm       =      sttm;
    output.nsim       =      nsimcomp;
    output.SimCount   =      SimCount;
    output.ValidCount =      ValidCount;
    output.jobindex   =      njobindex;
    output.tmelapsed  =      tmelapsed;
    output.lencat     =      lencat(1:nsimcomp,:);
    output.validcat   =      validcat;
    output.Mt         =      Mt;
    output.pixlev     =      Sresol;
    output.tresol     =      tresol;
    output.runtm      =      toc;
    %% save simulations if necessary
    %sttmstr = datestr(sttm,30);
    sttmstr = num2str(sttm);
    varnam = strcat(simfoldnam,'Sim_sttm',sttmstr,'_Maux',num2str(Maux),...
        '_Mpri',num2str(Mpri),'_KDEBKG',num2str(KDE_BKG),'_OmTyp',num2str(model.OmTyp),...
        '_btyp',num2str(btyp),'_isotyp',num2str(isotyp),'_srclen',num2str(srclen),'_mthres',num2str(mthres),...
        '_prodtyp',num2str(prodtyp),'_seqSpec',num2str(seqSpec),'_minFamSz',num2str(minFamSz),'_nprod',num2str(nprod),...
        '_tresol',num2str(tresol),'_jobindex',num2str(njobindex),'.mat');
    
    save(varnam,'output')
    
end
