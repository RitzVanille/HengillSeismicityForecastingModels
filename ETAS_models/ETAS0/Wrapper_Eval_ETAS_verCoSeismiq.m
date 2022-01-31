function [sim,modeleval] = Wrapper_Eval_ETAS_verCoSeismiq(SimUpdateDur,jobid)
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
% [sim,modeleval] = Wrapper_Eval_ETAS_verCoSeismiq(30,13)

% some parameters of obtaining smooth empirical pdf
%maxpossnum = 1000;
div        = 10;
combsim = 0;
comptyp = computer

baseoutfolder = './outputvars/';

%%
eqArea = 0;
savindvlik = 1;
rats_probs = 0;

catnam   =  'Coseismiq';
Regionnam = 'Iceland-Hengill-Active_pval0.1';
Prefix = strcat(catnam,'_',Regionnam);

cat0   =  importdata(strcat(Prefix,'.mat'));



simfoldnam = strcat(baseoutfolder,Prefix,'/SimComb_updateDur',num2str(SimUpdateDur),'/');
mkdir(simfoldnam)
evalfoldnam = strcat(baseoutfolder,Prefix,'/Eval_updateDur',num2str(SimUpdateDur),'/');
mkdir(evalfoldnam)

simpath  =  strcat(baseoutfolder,Prefix,'/Sim_updateDur',num2str(SimUpdateDur),'/');


tic

umodcomb = importdata(['./EvalList_',Regionnam,'.mat']);

curjobchunk = jobid;
%% import one simulation
%modnam      =    modlist(modid).name
for jobs = 1:length(curjobchunk)
    modid = curjobchunk(jobs);
    sttm        =    umodcomb.sttm(modid);
    sttmstr     = num2str(sttm);
    Maux        =    umodcomb.Maux(modid);
    Mpri        =    umodcomb.Mpri(modid);
    OmTyp       =    umodcomb.OmTyp(modid);
    btyp        =    umodcomb.btyp(modid); % this is temporary
    KDE_BKG     =    umodcomb.KDE_BKG(modid);
    isotyp      =    umodcomb.isotyp(modid);
    tresol      =    umodcomb.tresol(modid);
    srclen      =    umodcomb.srclen(modid);
    nprod      =    umodcomb.nprod(modid);
    prodtyp      =    umodcomb.prodtyp(modid);
    mthres      =    umodcomb.mthres(modid);
    minFamSz = umodcomb.minFamSz(modid);
    seqSpec = umodcomb.seqSpec(modid);
    M0          =    Mpri;
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
    end
    %%
    switch prodtyp
        case 0
            ProdText  =  'Parameteric' % temporary
        case 1
            ProdText = 'Bilinear Parametric with same spatial kernel'
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
    
    varnam = strcat(simpath,'Sim_sttm',sttmstr,'_Maux',num2str(Maux),...
        '_Mpri',num2str(Mpri),'_KDEBKG',num2str(KDE_BKG),'_OmTyp',num2str(OmTyp),...
        '_btyp',num2str(btyp),'_isotyp',num2str(isotyp),'_srclen',num2str(srclen),'_mthres',num2str(mthres),...
        '_prodtyp',num2str(prodtyp),'_seqSpec',num2str(seqSpec),'_minFamSz',num2str(minFamSz),'_nprod',num2str(nprod),...
        '_tresol',num2str(tresol),'_jobindex*.mat');
    
    simlist     =  dir(varnam);
    
    
    %% store
    sim       =    importdata(strcat(simpath,simlist(1).name));
    tresol    =    sim.tresol;
    Sresol    =    sim.pixlev;
    Mt        =    sim.Mt;
    lencat    =    [];
    nsim      =    0;
    SimCount = cell([length(Sresol),length(Mt)]);
    
    check = 1;
    for l = 1:size(simlist,1)
        
        try
            sim       =    importdata(strcat(simpath,simlist(l).name));
            
            for i = 1:length(Sresol)
                for j = 1:length(Mt)
                    SimCount{i,j} = [SimCount{i,j}; sim.SimCount{i,j} , ...
                        repmat(l,[length(sim.SimCount{i,j}(:,1)),1])];
                end
            end
            
            
            lencat = [lencat;sim.lencat];
            nsim   = nsim + sim.nsim;
        catch
            check = 1;
        end
    end
    
    check = 1;
    
    %% sort and add the counts
    
    
    nSimCount = cell([length(Sresol),length(Mt)]);
    for k = 1:length(Sresol)
        for j = 1:length(Mt)
            SimCount{k,j} = sortrows(SimCount{k,j}, [3,4,1,5]);
            [C, iaf] = unique(SimCount{k,j}(:,[3,4,1]),'rows');
            [C, ial] = unique(SimCount{k,j}(:,[3,4,1]),'last','rows');
            ntot     = zeros(size(iaf));
            for i = 1:length(iaf)
                ntot(i) = sum(SimCount{k,j}(iaf(i):ial(i),2));
            end
            
            nSimCount{k,j} = [C,ntot];
            
            %if rats_probs == 1
            % this is to compute mean rate to make
            % plots
            [ugrid_rate{k,j},iaf]  = unique(nSimCount{k,j}(:,[1,2]),'rows');
            [ugrid_rate{k,j},ial]  = unique(nSimCount{k,j}(:,[1,2]),'rows','last');
            rat     = zeros(size(iaf));
            for i = 1:length(iaf)
                rat(i) = sum(nSimCount{k,j}(iaf(i):ial(i),4).*nSimCount{k,j}(iaf(i):ial(i),3))/nsim;
            end
            
            [LonGMat,LatGMat,areamat,grid_vec_format] = GlobRectGrid_ver3([Sresol(k),Sresol(k)],...
                [min(cat0.Cregion(:,1)),max(cat0.Cregion(:,1))],[min(cat0.Cregion(:,2)),max(cat0.Cregion(:,2))]);
            tformgrid = floor([grid_vec_format(:,1)/Sresol(k), grid_vec_format(:,2)/Sresol(k)]) + 1;
            [lia,locb] = ismember(tformgrid,ugrid_rate{k,j},'rows');
            ratgrid = nan(size(lia));
            ratgrid(lia) = rat(locb(lia));
            
            ratmat = vec2mat(ratgrid,size(LonGMat,1))';
            
            sim.LonGMat{k} = LonGMat;
            sim.LatGMat{k} = LatGMat;
            sim.areamat{k} = areamat;
            sim.ratmat{k,j} = ratmat;
            sim.tformgrid{k} = tformgrid;
            %end
            
            
        end
    end
    
    %% save the combined models
    
    sim.SimCount = nSimCount;
    sim.nsim     = nsim;
    sim.lencat   = lencat;
    try
        sim = rmfield(sim,'Allcat');
    catch
    end
    %sim.pixlev   = [1:8]; % can be removed later
    varnam0 = strcat(simfoldnam,'Sim_sttm',sttmstr,'_Maux',num2str(Maux),...
        '_Mpri',num2str(Mpri),'_KDEBKG',num2str(KDE_BKG),'_OmTyp',num2str(OmTyp),...
        '_btyp',num2str(btyp),'_isotyp',num2str(isotyp),'_srclen',num2str(srclen),'_mthres',num2str(mthres),...
        '_prodtyp',num2str(prodtyp),'_seqSpec',num2str(seqSpec),'_minFamSz',num2str(minFamSz),'_nprod',num2str(nprod),...
        '_tresol',num2str(tresol),'.mat');
    
    save(varnam0,'sim');
    
    
    %% Evaluate Models
    wlev = 1/(sim.nsim+1);
    disp('-----------------------------------------------------------------')
    disp(jobs)
    disp(datestr(sim.sttm))
    disp('-----------------------------------------------------------------')
    %[loglik,logliksum] = EvaluateModels_ver3(sim,wlev);
    
    [loglik,logliksum] = EvaluateModels_ver601_verdeg(sim,div);
    
    
    if savindvlik == 1
        modeleval.loglik    =       loglik;
    end
    modeleval.logliksum  =       logliksum;
    modeleval.pars       =     sim.pars;
    modeleval.Mt         =      sim.Mt;
    modeleval.pixlev     =      sim.pixlev;
    modeleval.nsim       =      sim.nsim;
    modeleval.sttm       =      sim.sttm;
    modeleval.div        =      div;
    modeleval.wlev       =      wlev;
    %modeleval.tmelapsed  =      totrumtm;
    modeleval.validcat   =      sim.validcat;
    modeleval.ValidCount   =    sim.ValidCount;
    
    varnam = strcat(evalfoldnam,'Eval_sttm',sttmstr,'_Maux',num2str(Maux),...
        '_Mpri',num2str(Mpri),'_KDEBKG',num2str(KDE_BKG),'_OmTyp',num2str(OmTyp),...
        '_btyp',num2str(btyp),'_isotyp',num2str(isotyp),'_srclen',num2str(srclen),'_mthres',num2str(mthres),...
        '_prodtyp',num2str(prodtyp),'_seqSpec',num2str(seqSpec),'_minFamSz',num2str(minFamSz),'_nprod',num2str(nprod),...
        '_tresol',num2str(tresol),'.mat');
    
    save(varnam,'modeleval');
    toc
end


end
