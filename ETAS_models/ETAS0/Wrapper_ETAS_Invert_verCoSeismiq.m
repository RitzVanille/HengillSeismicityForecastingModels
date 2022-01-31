function [output] = Wrapper_ETAS_Invert_verCoSeismiq(maxtm,Maux,Mpri,KDE_BKG,OmTyp,btyp,isotyp,prodtyp,...
    srclen,dur,mthres,seqSpec,minFamSize,finiteInteg,id)
%% citations

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


%}

%%

%{
   
    
    % for iceland
    maxtm = inf;
    Maux = 0.3;
    Mpri = 0.3;
    KDE_BKG = 1; % can be 0
    OmTyp = 1;
    btyp  = 1;
    isotyp = 1;
    prodtyp = 0;
    srclen = 50;
    dur = 30; % duration of the testing period.
    mthres = 0;
    finiteInteg = 0;
    seqSpec = 0;
    minFamSize = 0;
    id = 381; % id of the testing period. min value for id is 1 and max
    %value is N, where N is the maximum number of testing period. 
    
    
    [output] = Wrapper_ETAS_Invert_verCoSeismiq(maxtm,Maux,Mpri,KDE_BKG,...
     OmTyp,btyp,isotyp,prodtyp,srclen,dur,mthres,seqSpec,minFamSize,finiteInteg,id)


    %}
    
    %%
    tmbegin = tic;
    rng(id)
    
    baseoutfolder = './outputvars/';
    %baseoutfolder = '/media/shyam/Seagate Backup Plus Drive/Shyam/datadump/KDE_BKG/outputvars/';
    %baseoutfolder = 'G:/Shyam/datadump/KDE_BKG/outputvars/';
    %baseoutfolder = 'T:/Work/datadump/KDE_BKG/outputvars/';
    %baseoutfolder = 'R:/datadump/KDE_BKG/outputvars/';
    %% parameters that remain fixed in this experiment
    
    %minK     = 45; % can be set by trial and error
    SMkernel = 1;
    anisodur = 365;
    mintrigprob = 0.05;
    
    usecentroid = 0;
    savtrigprob = 0;
    curcat = 0;
    minK = 2;
    maxK     = minK + 100;
    maxdist  = 1000; % number of nn for each earthquake
    maxglobvars = 30*10^6;
    nn       = 1;
    M0       = Mpri;
    initprodK = 0.27;
    initproda = 2;
    tmpars = [tmbegin,maxtm];
    
    %% options for inversion
    globoptions.options = optimoptions('fmincon','Algorithm','interior-point','Display','off',...
        'TolFun',1e-2,'TolX',1e-4,'MaxFunEvals',100*50,'MaxIter',100);
    
    globoptions.tmpars = tmpars;
    globoptions.maxiter = 100;
    globoptions.tolpars = 10^-2;
    globoptions.tolLL = 0;
    globoptions.toltrigprob = 10;
    globoptions.numloop1 = 5;
    %%
    % get cat props from cat id
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
    %%
    
    cat0   =  importdata(strcat('./',Prefix,'.mat'));
    
    invfoldnam = strcat(baseoutfolder,Prefix,'/Inv_updateDur',num2str(dur),'/');
    mkdir(invfoldnam)
    
    styraux = cat0.styraux;
    
    
    
    %% init
        initpars = [];
        priormodel = [];
        sttm = [];
    %%
    dm = cat0.dm;
    mc = cat0.mc;
    if M0 < mc
        Mpri = mc;
        M0 = mc;
    end
    if Maux < min(cat0.auxcat.magnitude)
        Maux = min(cat0.auxcat.magnitude);
    end
    
    Aearth = cat0.Aearth;
    
    styrpri = cat0.styrpri;
    if styrpri<2100
        styrpri  =  datenum(styrpri,1,1);
    else
        styrpri = cat0.styrpri;
    end
    
    
    pricat = [cat0.pricat.datenum,cat0.pricat.longitude,cat0.pricat.latitude,cat0.pricat.magnitude];
    [pricat(:,4)]  =  BinMags(pricat(:,4),0,dm);
    auxcat = [cat0.auxcat.datenum,cat0.auxcat.longitude,cat0.auxcat.latitude,cat0.auxcat.magnitude];
    [auxcat(:,4)]  =  BinMags(auxcat(:,4),0,dm);
    
    
    %%
    if KDE_BKG == 1
        SMkernel = 1;
    elseif KDE_BKG == 2
        SMkernel = 4;
    elseif KDE_BKG == 3
        SMkernel = 5;
    end
    
    
    %% define the starting times of the testing periods
    if isempty(sttm)
        endsttm = max(cat0.pricat.datenum);
        sttmlist = [begsttm:dur:endsttm];
        if id > length(sttmlist)
            disp('Setting id to maxid')
            id = length(sttmlist)
        end
        sttm     = sttmlist(id);
        smHyper = [maxdist,maxK,cat0.mind,SMkernel,minK];
    else
        if KDE_BKG~=0
            smHyper = [priormodel.maxdist,priormodel.SMpars(1)+6,...
                cat0.mind,priormodel.SMkernel,priormodel.SMpars(1)-4];
            SMkernel = priormodel.SMkernel;
        else
            smHyper = [];
        end
    end
    
    
    if size(curcat,1)>1
        sttm = max(auxcat(:,1));
    end
    
    disp(Prefix)
    disp(datestr(sttm))
    
    %% prepare the catalog
    
    catsrc =  auxcat(auxcat(:,1) >= styraux & auxcat(:,1) < sttm & auxcat(:,4)>= Maux,:);
    catsrc_special =  auxcat(auxcat(:,1) >= styraux & auxcat(:,1) < sttm & auxcat(:,4)>= Mpri,:);
    if  finiteInteg == 1 || ~isempty(cat0.Cregion)
        [IN] = inpolygon(catsrc(:,2),catsrc(:,3),cat0.Cregion(:,1),cat0.Cregion(:,2));
        catsrc = catsrc(IN,:);
    end
    
    catsrctab = cat0.auxcat(auxcat(:,1) >= styraux & auxcat(:,1) < sttm & auxcat(:,4)>= Maux,:);
    catrec =  pricat(pricat(:,1) >= styrpri & pricat(:,1) < sttm & pricat(:,4)>= Mpri,:);
    catrectab = cat0.pricat(pricat(:,1) >= styrpri & pricat(:,1) < sttm & pricat(:,4)>= Mpri,:);
    
    timewin     =  [min(catrec(:,1)),max(catsrc(:,1))];
    indcatrec =  (ismember(catsrc,catrec,'rows'));
    indcatrec_special =  find(ismember(catsrc_special,catrec,'rows'));
    %% change the maxdist for the SMkernel based on the size of the catalog
    maxdist = min(round(maxglobvars/length(catrec)),maxdist);
    if maxdist > length(catrec(:,1))
        maxdist = length(catrec(:,1));
    end
    smHyper(1) =  maxdist;
    %% b-value
    [beta]   =  Tinti(catrec(:,4),min(catrec(:,4)),dm);
    
    
    %% Model description
    
    switch seqSpec
        case 0
            SEQ_Text = 'Global Parameters' % space and temporally homogenous poisson process
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
        case 3
            KDE_Text = 'SVTV_PP' % space varying temporally varying poisson process
    end
    %%
    switch prodtyp
        case 0
            ProdText  =  'Parameteric' % temporary
        case 1
            ProdText = 'Bilinear Parametric with same spatial kernel'
            mthres
        case 2
            ProdText = 'Bilinear Parametric with different spatial kernel'
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
        case 201
            OmText = 'Omori + Exponential'
        case 202
            OmText = 'quadratic p(m)'
            
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
            %{
        case 2
            isoText = 'Anisotropic' % temporary
            % load the anistropic clusters
            %AnisoClusters = importdata('../inpvars/AnisoClust_ANSS_Globe.mat');
            anisoInpcat = cat0;
            anisoInpcat.auxcat(anisoInpcat.auxcat.datenum >= sttm,:) = [];
            anisoInpcat.pricat(anisoInpcat.pricat.datenum >= sttm,:) = [];
            disp('Building Anisotropic Cluster...')
            [AnisoClusters] = Wrapper_NN_Aniso_SKern_ithEQ_ver5(anisoInpcat,anisodur,mthres,[]);
            AnisoClusters.largeshock = anisoInpcat.pricat(anisoInpcat.pricat.magnitude >= mthres , :);
            largeshocks = [AnisoClusters.largeshock.datenum,AnisoClusters.largeshock.longitude,...
                AnisoClusters.largeshock.latitude,AnisoClusters.largeshock.magnitude];
            % find the index of the largeshocks in the catsrctab
            [Lia,Locb] = ismember(catsrc,largeshocks,'rows');
            AftClust = AnisoClusters.After_Cluster(Locb(Lia));
            PriorClust = AnisoClusters.Prior_Cluster(Locb(Lia));
            
            LargeMainId = find(Lia);
            largeshocks = largeshocks(Locb(Lia),:);
        case 3
            isoText = 'Isotropic 2.0' % temporary
        case 4
            isoText = 'Magnitude dependent rho' % temporary
            %}
    end
    
    %datestr(sttm)
    %% compute the globvars
    %{
    if isotyp == 2
        disp('computing Anisotropic globvars...')
        [anisoglobvars.globvars_After] = GlobVars_LargeMain(AftClust,catrec,largeshocks,LargeMainId,srclen);
        [anisoglobvars.globvars_Prior] = GlobVars_LargeMain(PriorClust,catrec,largeshocks,LargeMainId,srclen);
        %anisoglobvars.globvars_After = globvars_After_LM;
        %anisoglobvars.globvars_Prior = globvars_prior_LM;
        anisoglobvars.largeshocks = largeshocks;
        anisoglobvars.LargeMainId = LargeMainId;
    end
    %}
    
    
        disp('Computing Globvars for EM...')
       
            [globvars]  =  GlobVars(catsrc,catrec,srclen,0);
        
        
        % globvars for the anisotropic clusters
        
        globvars(globvars(:,1) == 0,1) = 10^-2;
        
            disp('Computing globvars for smoothing...')
            [globvars_bkg,maxrad] = GlobvarsDistSSM_ver3(catrec(:,2:3),maxdist);
        
    
    
    %% initial parameters
    prodTable = [];
    lhat_corr = [];
    ProdPars  = [];
    SMpars = [];
    Adaptive_BW = [];
    mu_svar = [];
    
    %% inversion
    if KDE_BKG == 0 && OmTyp == 1 && btyp == 1 && isotyp == 1 && prodtyp == 0 && seqSpec == 0
        % null model
        mthres = 0;
        [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,Invcomp] =  Invert_EtasIso_Omori1(catsrc,M0,globvars,...
            Aearth,timewin,catrec,initpars,globoptions);
        partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
            'gamma','c0','p0','tau0'});
        %output.trigprob = trigprob(trigprob(:,1)>=mintrigprob,:);
    elseif  KDE_BKG == 1 && OmTyp == 1 && btyp == 1 && isotyp == 1 && prodtyp == 0 && seqSpec == 0
        % very promising 12 august 2021: current richterx model
        mthres = 0;
        [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp] =  ...
            Invert_EtasIso_Omori1_KDEBKG_ver2(catsrc,M0,globvars,Aearth,timewin,catrec,initpars,priormodel,smHyper,...
            globoptions,globvars_bkg,maxrad);
        partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
            'gamma','c0','p0','tau0'});
      %{  
    elseif  KDE_BKG == 1 && OmTyp == 2 && btyp == 1 && isotyp == 1 && prodtyp == 0 && seqSpec == 0
        % not so promising
        check = 1;
        [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
            Invert_EtasIso_Omori2_KDEBKG(catsrc,M0,globvars,Aearth,timewin,catrec,initpars,priormodel,smHyper,...
            globoptions,globvars_bkg,maxrad);
        partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
            'gamma','c0','c1','p0','p1','tau0'});
        
    elseif  KDE_BKG == 1 && OmTyp == 2 && btyp == 2 && isotyp == 1 && prodtyp == 0 && seqSpec == 0
        %initpars = [-8.67,0.05,3.38,2.04,1.21,1.34,-4.74,0.54,1.05,0.01,3.37,2.32,1.75,2.98];
        check = 1;
        mthres = 0;
        [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
            Invert_EtasIso_Omori2_b2_KDEBKG_ver3(catsrc,M0,globvars,Aearth,timewin,catrec,...
            initpars,priormodel,smHyper,globoptions,globvars_bkg,maxrad);
        partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
            'gamma','c0','c1','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2'});
        
    elseif  KDE_BKG == 1 && OmTyp == 2 && btyp == 2 && isotyp == 1 && prodtyp == 1 && seqSpec == 0
        
        [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
            Invert_EtasIso_Omori2_b2_KDEBKG_Prod1(catsrc,M0,globvars,Aearth,timewin,catrec,...
            initpars,priormodel,smHyper,globoptions,globvars_bkg,maxrad,mthres);
        partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
            'gamma','c0','c1','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2','a2_unNorm'});
        partable.K2Norm = ProdPars(3);
        partable.a2Norm = ProdPars(4);
        
    elseif  KDE_BKG == 1 && OmTyp == 2 && btyp == 2 && isotyp == 1 && prodtyp == 3 && seqSpec == 0
        
        [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
            Invert_EtasIso_Omori2_b2_KDEBKG_Prod3(catsrc,M0,globvars,Aearth,timewin,catrec,...
            initpars,priormodel,smHyper,globoptions,globvars_bkg,maxrad,mthres);
        partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
            'gamma','c0','c1','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2','a2_unNorm'});
        partable.a2Norm = ProdPars(3);
        
    elseif  KDE_BKG == 1 && OmTyp == 2 && btyp == 2 && isotyp == 1 && prodtyp == 301 && seqSpec == 0
        
        [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
            Invert_EtasIso_Omori2_b2_KDEBKG_Prod301(catsrc,M0,globvars,Aearth,timewin,catrec,...
            initpars,priormodel,smHyper,globoptions,globvars_bkg,maxrad);
        partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
            'gamma','c0','c1','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2','a2_unNorm'});
        partable.a2Norm = ProdPars(3);
        partable.Mcritical = ProdPars(4);
        partable.acritical = ProdPars(5);
        
        
    elseif  KDE_BKG == 1 && OmTyp == 2 && btyp == 2 && isotyp == 1 && prodtyp == 303 && seqSpec == 0
        
        [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
            Invert_EtasIso_Omori2_b2_KDEBKG_Prod303(catsrc,M0,globvars,Aearth,timewin,catrec,...
            initpars,priormodel,smHyper,globoptions,globvars_bkg,maxrad);
        partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
            'gamma','c0','c1','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2','a2_unNorm','a3_unNorm'});
        partable.a2Norm = ProdPars(3);
        partable.a3Norm = ProdPars(4);
        
    elseif  KDE_BKG == 1 && OmTyp == 2 && btyp == 2 && isotyp == 1 && prodtyp == 304 && seqSpec == 0
        
        [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
            Invert_EtasIso_Omori2_b2_KDEBKG_Prod304(catsrc,M0,globvars,Aearth,timewin,catrec,...
            initpars,priormodel,smHyper,globoptions,globvars_bkg,maxrad);
        partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
            'gamma','c0','c1','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2','a2_unNorm'});
        partable.a2Norm = ProdPars(3);
        
        
        %partable.a3Norm = ProdPars(4);
        
    elseif  KDE_BKG == 3 && OmTyp == 2 && btyp == 2 && isotyp == 1 && prodtyp == 301 && seqSpec == 0
        % all inversions need to be finished and need to be tested
        %maxrad = [];
        %         initpars = [-6.1067   -0.3416   -0.2883   -0.5711    0.6799    1.1294   -4.8079...
        %             0.5248 0.3376    0.1778    3.2933    2.5813    1.8178    3.4202...
        %             0.2437];
        globvars_pri = globvars(globvars(:,4) >= min(indcatrec_special),:);
        globvars_pri(:,4) = globvars_pri(:,4) - min(indcatrec_special) + 1;
        [pars,LLall,trigprob,lhat_estep,IP,IP_aux,AIP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
            Invert_EtasIso_Omori2_b2_KDEBKG3_Prod301(catsrc,M0,globvars_pri,Aearth,timewin,catrec,...
            initpars,priormodel,smHyper,globoptions,globvars_bkg,globvars_bkg_tdep,maxrad,globvars,catsrc_special);
        partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
            'gamma','c0','c1','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2','a2_unNorm'});
        partable.a2Norm = ProdPars(3);
        partable.Mcritical = ProdPars(4);
        partable.acritical = ProdPars(5);
        output.AIP = AIP;
        output.IP_aux = IP_aux;
        
        elseif  KDE_BKG == 1 && OmTyp == 2 && btyp == 2 && isotyp == 4 && prodtyp == 3 && seqSpec == 0
        
        [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
            Invert_Etas_Iso4_Omori2_b2_KDEBKG_Prod3(catsrc,M0,globvars,Aearth,timewin,catrec,...
            initpars,priormodel,smHyper,globoptions,globvars_bkg,maxrad);
        
        partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,...
            'gamma','c0','c1','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2','a2_unNorm','q0','q1'});
        partable.a2Norm = ProdPars(3);
        
        
        
        elseif  KDE_BKG == 1 && OmTyp == 2 && btyp == 2 && isotyp == 1 && prodtyp == 302 && seqSpec == 0
        
        [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
            Invert_EtasIso_Omori2_b2_KDEBKG_Prod302(catsrc,M0,globvars,Aearth,timewin,catrec,...
            initpars,priormodel,smHyper,globoptions,globvars_bkg,maxrad);
        partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
            'gamma','c0','c1','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2','a2_unNorm','M0_min','Del'});
        partable.a2Norm = ProdPars(3);
        partable.Mcritical = ProdPars(6);
        partable.acritical = ProdPars(7);
        
        
        
        elseif  KDE_BKG == 1 && OmTyp == 202 && btyp == 2 && isotyp == 1 && prodtyp == 301 && seqSpec == 0
        
        [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
            Invert_EtasIso_Omori202_b2_KDEBKG_Prod301(catsrc,M0,globvars,Aearth,timewin,catrec,...
            initpars,priormodel,smHyper,globoptions,globvars_bkg,maxrad);
        partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
            'gamma','c0','c1','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2','a2_unNorm','p2'});
        partable.a2Norm = ProdPars(3);
        partable.Mcritical = ProdPars(4);
        partable.acritical = ProdPars(5);
        
        
             elseif  KDE_BKG == 1 && OmTyp == 201 && btyp == 2 && isotyp == 1 && prodtyp == 301 && seqSpec == 0
        
                [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
                    Invert_EtasIso_Omori201_b2_KDEBKG_Prod301(catsrc,M0,globvars,Aearth,timewin,catrec,...
                    initpars,priormodel,smHyper,globoptions,globvars_bkg,maxrad);
                partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
                    'gamma','c0','c1','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2','a2_unNorm'});
                partable.a2Norm = ProdPars(3);
                partable.Mcritical = ProdPars(4);
                partable.acritical = ProdPars(5);
        
            elseif  KDE_BKG == 1 && OmTyp == 1 && btyp == 1 && isotyp == 1 && prodtyp == 0 && seqSpec == 1
                % we need a prior model that corresponds to the all config same
                % except sequence specificity
                if isstruct(init)
                    priormodel = init;
                elseif init == 1
                    disp('Prior Model already assigned')
                else
                    priormodel = importdata(strcat(invfoldnam,'Pars_sttm',num2str(sttm),'_Maux',num2str(Maux),...
                        '_Mpri',num2str(Mpri),'_KDEBKG',num2str(KDE_BKG),'_OmTyp',num2str(OmTyp),'_btyp',num2str(btyp)...
                        ,'_isotyp',num2str(isotyp),'_srclen',num2str(srclen),'_mthres',num2str(mthres),'_prodtyp',num2str(prodtyp),...
                        '_seqSpec',num2str(0),'_minFamSz',num2str(0),'.mat'));
                end
                globoptions.tolpars = 10^-1.5;
                globoptions.maxiter = 30;
                globoptions.maxiterout = 100;
        
                priormodel.Aearth  = Aearth;
                priormodel.timewin  = timewin;
                priormodel.smHyper  = smHyper;
                priormodel.maxrad  = maxrad;
                mthres = 0;
                [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,EqFamilies,FamilySz,parevol,Invcomp] =  ...
                    Invert_EtasIso_Omori1_KDEBKG_seqSpec_outer(catsrc,catrec,globvars,priormodel,globoptions,globvars_bkg,minFamSize);
                partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
                    'gamma','c0','p0','tau0'});
                output.parevol = parevol;
        
            elseif  KDE_BKG == 1 && OmTyp == 2 && btyp == 2 && isotyp == 1 && prodtyp == 0 && seqSpec == 1
                % we need a prior model that corresponds to the all config same
                % except sequence specificity
                if isstruct(init)
                    priormodel = init;
                elseif init == 1
                    disp('Prior Model already assigned')
                else
                    priormodel = importdata(strcat(invfoldnam,'Pars_sttm',num2str(sttm),'_Maux',num2str(Maux),...
                        '_Mpri',num2str(Mpri),'_KDEBKG',num2str(KDE_BKG),'_OmTyp',num2str(OmTyp),'_btyp',num2str(btyp)...
                        ,'_isotyp',num2str(isotyp),'_srclen',num2str(srclen),'_mthres',num2str(mthres),'_prodtyp',num2str(prodtyp),...
                        '_seqSpec',num2str(0),'_minFamSz',num2str(0),'.mat'));
                end
                globoptions.tolpars = 10^-1.5;
                globoptions.maxiter = 30;
                globoptions.maxiterout = 100;
        
                priormodel.Aearth  = Aearth;
                priormodel.timewin  = timewin;
                priormodel.smHyper  = smHyper;
                priormodel.maxrad  = maxrad;
                priormodel.indcatrec = indcatrec;
                mthres = 0;
        
                [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,...
                    Adaptive_BW,mu_svar,EqFamilies,FamilySz,parevol,Invcomp] =  ...
                    ...
                    Invert_EtasIso_Omori2_b2_KDEBKG_seqSpec_outer(catsrc,catrec,globvars,...
                    priormodel,globoptions,globvars_bkg,minFamSize);
        
                partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
                    'gamma','c0','p0','tau0'});
                output.parevol = parevol;
            elseif  KDE_BKG == 1 && OmTyp == 1 && btyp == 1 && isotyp == 1 && prodtyp == 1
                % very promising
                %mthres = 7;
                [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp] =  ...
                    Invert_EtasIso_Omori1_KDEBKG_Prod1(catsrc,M0,globvars,Aearth,timewin,catrec,initpars,priormodel,smHyper,...
                    globoptions,globvars_bkg,maxrad,mthres);
                partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
                    'gamma','c0','p0','tau0','a2_unNorm'});
                partable.K2Norm = ProdPars(3);
                partable.a2Norm = ProdPars(4);
        
            elseif  KDE_BKG == 1 && OmTyp == 1 && btyp == 2 && isotyp == 1 && prodtyp == 1
                % very promising
                %mthres = 7;
                [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp] =  ...
                    Invert_EtasIso_Omori1_KDEBKG_Prod1_b2(catsrc,M0,globvars,Aearth,timewin,catrec,initpars,priormodel,smHyper,...
                    globoptions,globvars_bkg,maxrad,mthres);
                partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
                    'gamma','c0','p0','tau0','a2_unNorm','beta_bkg','beta_aft1','beta_aft2'});
                partable.K2Norm = ProdPars(3);
                partable.a2Norm = ProdPars(4);
        
            elseif  KDE_BKG == 1 && OmTyp == 1 && btyp == 1 && isotyp == 1 && prodtyp == 2
                % very promising
                %mthres = 7;
                [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp] =  ...
                    Invert_EtasIso_Omori1_KDEBKG_Prod2(catsrc,M0,globvars,Aearth,timewin,catrec,initpars,priormodel,smHyper,...
                    globoptions,globvars_bkg,maxrad,mthres);
                partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
                    'gamma','c0','p0','tau0','a2_unNorm','d1','rho1','gamma1'});
                partable.K2Norm = ProdPars(3);
                partable.a2Norm = ProdPars(4);
        
        output.trigprob = trigprob(trigprob(:,1)>=mintrigprob,:);
        elseif  KDE_BKG == 1 && OmTyp == 1 && btyp == 2 && isotyp == 1
            % tbi
            [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp] =  ...
                Invert_EtasIso_Omori1_b2_KDEBKG_corr(catsrc,M0,globvars,Aearth,timewin,catrec,initpars,priormodel,smHyper,...
                globoptions,globvars_bkg,maxrad);
            partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
                'gamma','c0','p0','tau0','beta_bkg','beta_aft1','beta_aft2'});
        output.trigprob = trigprob(trigprob(:,1)>=mintrigprob,:);
            %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif  KDE_BKG == 0 && OmTyp == 2 && btyp == 1
        
            [pars,LLall,~,lhat_estep,IP,prodTable,lhat_corr,ProdPars,Invcomp]   =  ...
                Invert_EtasIso_Omori2(catsrc,M0,globvars,Aearth,timewin,catrec,initpars,tmpars);
            partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
                'gamma','c0','c1','p0','p1','tau0'});
        
        elseif  KDE_BKG == 1 && OmTyp == 2 && btyp == 1 && isotyp == 1
            % not so promising
            check = 1;
            [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
                Invert_EtasIso_Omori2_KDEBKG(catsrc,M0,globvars,Aearth,timewin,catrec,initpars,priormodel,smHyper,...
                globoptions,globvars_bkg,maxrad);
            partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
                'gamma','c0','c1','p0','p1','tau0'});
        output.trigprob = trigprob(trigprob(:,1)>=mintrigprob,:);
        
        
        
        
            elseif  KDE_BKG == 2 && OmTyp == 2 && btyp == 2 && isotyp == 1 && prodtyp == 0
                %initpars = [-8.3631   -0.3345    1.9939    2.0624    0.9385    0.7404   -4.0932    0.3029 0.1720    0.1406    3.6932    2.2646    2.0059    3.2007];
                check = 1;
                mthres = 0;
                [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
                    Invert_EtasIso_Omori2_b2_KDEBKG2(catsrc,M0,globvars,Aearth,timewin,catrec,...
                    initpars,priormodel,smHyper,globoptions,globvars_bkg,maxrad);
                partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
                    'gamma','c0','c1','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2'});
        
            elseif  KDE_BKG == 2 && OmTyp == 2 && btyp == 2 && isotyp == 1 && prodtyp == 1
                %initpars = [-8.3631   -0.3345    1.9939    2.0624    0.9385    0.7404   -4.0932    0.3029 0.1720    0.1406    3.6932    2.2646    2.0059    3.2007 2.7];
        
                [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
                    Invert_EtasIso_Omori2_b2_KDEBKG2_Prod1(catsrc,M0,globvars,Aearth,timewin,catrec,...
                    initpars,priormodel,smHyper,globoptions,globvars_bkg,maxrad,mthres);
                partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
                    'gamma','c0','c1','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2','a2_unNorm'});
                partable.K2Norm = ProdPars(3);
                partable.a2Norm = ProdPars(4);
        output.trigprob = trigprob(trigprob(:,1)>=mintrigprob,:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif  KDE_BKG == 0 && OmTyp == 3 && btyp == 1
        
            [pars,LLall,~,lhat_estep,IP,prodTable,lhat_corr,ProdPars]   =  ...
                Invert_EtasIso_Omori3(cat,M0,globvars,Aearth,timewin,catrec,initpars);
            partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
                'gamma','c0','p0','p1','tau0'});
        
        elseif  KDE_BKG == 1 && OmTyp == 3 && btyp == 2 && isotyp == 1
        
            [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
                Invert_EtasIso_Omori3_b2_KDEBKG(catsrc,M0,globvars,Aearth,timewin,catrec,initpars,priormodel,smHyper,globoptions,globvars_bkg,maxrad);
            partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
                'gamma','c0','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2'});
        output.trigprob = trigprob(trigprob(:,1)>=mintrigprob,:);
            elseif KDE_BKG == 1 && OmTyp == 1 && btyp == 1 && isotyp == 2 && prodtyp == 0
                %initpars = importdata('initpars.mat')
        
                [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
                    Invert_Etas_Iso2_Omori1_b1_KDEBKG(catsrc,M0,globvars,Aearth,timewin,catrec,initpars,priormodel,smHyper,...
                    globoptions,globvars_bkg,maxrad,anisoglobvars);
                partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
                    'gamma','c0','p0','tau0','d1','rho1','d2','rho2','w1','K2_unNorm','a2_unNorm'});
                partable.anisodur = anisodur;
        
                partable.K2Norm = ProdPars(3);
                partable.a2Norm = ProdPars(4);
                output.Anisocluster = AnisoClusters;
        
        
            elseif KDE_BKG == 1 && OmTyp == 2 && btyp == 2 && isotyp == 2 && prodtyp == 0
                %initpars = importdata('initpars.mat')
        
                [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
                    Invert_Etas_Iso2_Omori2_b2_KDEBKG_ver2(catsrc,M0,globvars,Aearth,timewin,catrec,initpars,priormodel,smHyper,...
                    globoptions,globvars_bkg,maxrad,anisoglobvars);
                partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
                    'gamma','c0','c1','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2',...
                    'd1','rho1','d2','rho2','w1','K2_unNorm','a2_unNorm'});
                partable.anisodur = anisodur;
                partable.Mthres   = mthres;
                partable.K2Norm = ProdPars(3);
                partable.a2Norm = ProdPars(4);
                output.Anisocluster = AnisoClusters;
        output.trigprob = trigprob(trigprob(:,1)>=mintrigprob,:);
        
        elseif  KDE_BKG == 1 && OmTyp == 2 && btyp == 2 && isotyp == 3
            %initpars = [-8.67,0.05,3.38,2.04,1.21,1.34,-4.74,0.54,1.05,0.01,3.37,2.32,1.75,2.98];
            check = 1;
            globvars(:,1) = sqrt(globvars(:,1));
            [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp]   =  ...
                Invert_Etas_Iso3_Omori2_b2_KDEBKG_corr(catsrc,M0,globvars,Aearth,timewin,catrec,...
                initpars,priormodel,smHyper,tmpars,globvars_bkg,maxrad);
            partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho',...
                'gamma','c0','c1','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2'});
            %output.trigprob = trigprob(trigprob(:,1)>=mintrigprob,:);
        
        elseif  KDE_BKG == 1 && OmTyp == 3 && btyp == 2
            [pars,LLall,~,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar]   =  ...
                Invert_EtasIso_Omori3_b2_KDEBKG(cat,M0,globvars,Aearth,timewin,catrec,initpars,priormodel);
            partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho','gamma',...
                'c0','p0','p1','tau0','beta_bkg','beta_aft1','beta_aft2'});
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % elseif  KDE_BKG == 0 && OmTyp == 4 && btyp == 1
            % elseif  KDE_BKG == 1 && OmTyp == 4 && btyp == 1
            %     % not so promising
            %     [pars,LLall,~,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar]   =  ...
            %         Invert_EtasIso_Omori4_KDEBKG(cat,M0,globvars,Aearth,timewin,catrec,initpars,priormodel);
            %     partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho','gamma','c0','c1','p0','p1','tau0','tau1'});
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif  KDE_BKG == 0 && OmTyp == 0 && btyp == 1
        
            [pars,LLall,~,lhat_estep,IP,prodTable,lhat_corr,ProdPars]   =  ...
                Invert_EtasIso_Omori0(cat,M0,globvars,Aearth,timewin,catrec,initpars);
            partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho','gamma','c0','p0'});
        
        elseif  KDE_BKG == 1 && OmTyp == 0 && btyp == 1
            % needs to be evaluated could be promising
            [pars,LLall,~,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar]   =  ...
                Invert_EtasIso_Omori0_KDEBKG(cat,M0,globvars,Aearth,timewin,catrec,initpars,priormodel);
            partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho','gamma','c0','p0'});
        
        elseif  KDE_BKG == 1 && OmTyp == 0 && btyp == 2
            % needs to be evaluated could be promising
            [pars,LLall,~,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar]   =  ...
                Invert_EtasIso_Omori0_b2_KDEBKG(cat,M0,globvars,Aearth,timewin,catrec,initpars,priormodel);
            partable = array2table(pars,'VariableNames',{'mu','K_unNorm','a_unNorm','d' ,'rho','gamma',...
                'c0','p0','beta_bkg','beta_aft1','beta_aft2'});
        
        %}
    end
    
    
    
    %% find the smoothing radius for finding 1 background earthquake in the circle around each
    % event
    disp('Computing nearest neighbour distance now')
    [smrad] = SmoothRad(catrec(:,2:3),nn,IP,cat0.mind,[]);
    %% find branching ratio
    partable.betaGlob = repmat(beta,[size(pars,1),1]);
    partable.KNorm = ProdPars(:,1);
    partable.aNorm = ProdPars(:,2);
    partable.Nbkg = repmat(sum(IP),[size(pars,1),1]);
    partable.Mthres   = repmat(mthres,[size(pars,1),1]);
    n = (size(catrec,1)-sum(IP))/size(catrec,1);
    n1 = mean(lhat_corr(:,2));
    n2 = n1;
    partable.n = repmat(n,[size(pars,1),1]);
    partable.n1 = repmat(n1,[size(pars,1),1]);
    partable.n2 = repmat(n2,[size(pars,1),1]);
    %partable.reg = reg;
    % save the outputs in the specified folder
    if KDE_BKG == 1
        partable.SMpars = repmat(SMpars,[size(pars,1),1]);
    end
    partable.isotyp = repmat(isotyp,[size(pars,1),1]);
    partable.isotext = repmat(isoText,[size(pars,1),1]);
    partable.seqSpec = repmat(seqSpec,[size(pars,1),1]);
    partable.SEQ_Text = repmat(SEQ_Text,[size(pars,1),1]);
    partable.minFamSize = repmat(minFamSize,[size(pars,1),1]);
    partable.btyp = repmat(btyp,[size(pars,1),1]);
    partable.OmTyp = repmat(OmTyp,[size(pars,1),1]);
    partable.sttm = repmat(sttm,[size(pars,1),1]);
    partable.KDE_BKG = repmat(KDE_BKG,[size(pars,1),1]);
    partable.prodtyp = repmat(prodtyp,[size(pars,1),1]);
    partable.Mmax = repmat(Mmax,[size(pars,1),1]);
    partable.Maux = repmat(Maux,[size(pars,1),1]);
    partable.Mpri = repmat(Mpri,[size(pars,1),1]);
    partable.mc = repmat(mc,[size(pars,1),1]);
    partable.dm = repmat(dm,[size(pars,1),1]);
    partable.srclen = repmat(srclen,[size(pars,1),1]);
    partable.SMkernel = repmat(SMkernel,[size(pars,1),1]);
    partable.maxdist = repmat(maxdist,[size(pars,1),1]);
    partable.styraux    =    repmat(styraux,[size(pars,1),1]);
    partable.styrpri    =    repmat(styrpri,[size(pars,1),1]);
    partable.nn = repmat(nn,[size(pars,1),1]);
    partable.LL = repmat(LLall,[size(pars,1),1]);
    partable.OmText = repmat(OmText,[size(pars,1),1]);
    partable.bText = repmat(bText,[size(pars,1),1]);
    partable.ProdText = repmat(ProdText,[size(pars,1),1]);
    partable.KDE_Text = repmat(KDE_Text,[size(pars,1),1]);
    partable.catnam = repmat(catnam,[size(pars,1),1]);
    partable.Regionnam = repmat(Regionnam,[size(pars,1),1]);
    
    if seqSpec == 1
        partable.FamilySz = FamilySz;
        output.EQfam = EqFamilies;
        output.FamilySz = FamilySz;
    end
    
    % to be removed
    if savtrigprob == 1
        output.pricat = catrec;
    elseif savtrigprob == 2
        output.catrectab = catrectab;
        output.catsrctab = catsrctab;
        output.trigprob = trigprob;
        % this is just temp
    end
    output.partable   =   partable;
    output.seqSpec = seqSpec;
    output.SEQ_Text = SEQ_Text;
    output.minFamSize = minFamSize;
    output.isotyp = isotyp;
    output.isotext = isoText;
    output.betaGlob = beta;
    output.ProdPars = ProdPars;
    output.n = [n,n1,n2];
    output.Nbkg = sum(IP);
    output.pars= pars;
    output.btyp = btyp;
    output.OmTyp = OmTyp;
    output.sttm = sttm;
    output.KDE_BKG = KDE_BKG;
    output.Mmax = Mmax;
    output.Maux = Maux;
    output.Mpri = Mpri;
    output.mc = mc;
    output.dm = dm;
    output.srclen = srclen;
    output.SMkernel = SMkernel;
    output.maxdist = maxdist;
    output.styraux    =    styraux;
    output.styrpri    =    styrpri;
    output.nn = nn;
    output.LL = LLall;
    output.OmText = OmText;
    output.bText = bText;
    output.ProdText = ProdText;
    output.KDE_Text = KDE_Text;
    output.Regionnam = Regionnam;
    output.catnam = catnam;
    output.Regionnam = Regionnam;
    output.InvComp    =    Invcomp;
    output.dur        =    dur;
    output.id         =    id;
    output.IP         =     IP;
    output.NumAft     =     lhat_estep;
    output.smrad      =    smrad;
    output.prodTable  =  prodTable;
    output.lhat_corr  =  lhat_corr;
    output.prodtyp     =  prodtyp;
    output.SMpars     =  SMpars;
    output.Adaptive_BW =  Adaptive_BW;
    output.mu_svar    =  mu_svar;
    
    output.runtum     =    toc(tmbegin);
    
    
    %%
    if Invcomp == 1
        
        varnam = strcat(invfoldnam,'Pars_sttm',num2str(sttm),'_Maux',num2str(Maux),...
            '_Mpri',num2str(M0),'_KDEBKG',num2str(KDE_BKG),'_OmTyp',num2str(OmTyp),'_btyp',num2str(btyp),...
            '_isotyp',num2str(isotyp),'_srclen',num2str(srclen),'_mthres',num2str(mthres),'_prodtyp',num2str(prodtyp),...
            '_seqSpec',num2str(seqSpec),'_minFamSz',num2str(minFamSize),'.mat');
        varnam2 = strcat(invfoldnam,'_Incomp_','Pars_sttm',num2str(sttm),'_Maux',num2str(Maux),...
            '_Mpri',num2str(M0),'_KDEBKG',num2str(KDE_BKG),'_OmTyp',num2str(OmTyp),'_btyp',num2str(btyp),...
            '_isotyp',num2str(isotyp),'_srclen',num2str(srclen),'_mthres',num2str(mthres),'_prodtyp',num2str(prodtyp),...
            '_seqSpec',num2str(seqSpec),'_minFamSz',num2str(minFamSize),'.mat');
        delete(varnam2)
    else
        varnam=strcat(invfoldnam,'_Incomp_','Pars_sttm',num2str(sttm),'_Maux',num2str(Maux),...
            '_Mpri',num2str(M0),'_KDEBKG',num2str(KDE_BKG),'_OmTyp',num2str(OmTyp),'_btyp',num2str(btyp),...
            '_isotyp',num2str(isotyp),'_srclen',num2str(srclen),'_mthres',num2str(mthres),'_prodtyp',num2str(prodtyp),...
            '_seqSpec',num2str(seqSpec),'_minFamSz',num2str(minFamSize),'.mat');
    end
    save(varnam,'output');
    toc(tmbegin)
    check=1;
end