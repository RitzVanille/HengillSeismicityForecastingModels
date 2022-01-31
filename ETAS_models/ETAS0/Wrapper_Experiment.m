 %{
Papers in which these codes were developed

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
%% work flow for ETAS modelling


catnam   =  'Coseismiq';
Regionnam = 'Iceland-Hengill-Active_pval0.1';
Prefix = strcat(catnam,'_',Regionnam);
cat0   =  importdata(strcat('./',Prefix,'.mat'));
%% First invert the parameters for different testing periods
% for iceland
maxtm = inf; % set this parameter to time assigned in seconds for a job on Euler
Maux = 0.3; % Mc
Mpri = 0.3; %Mc
KDE_BKG = 0; % can be 0 or 1, if 0 it leads to the basic ETAS model which does not have spatial variation of background rate, if 1 with spatially varying background rate
% some inputs that stay fixed in this project but are need in other
% projects.
%---------------------------------------------------
OmTyp = 1;
btyp  = 1;
isotyp = 1;
prodtyp = 0;
finiteInteg = 0;
seqSpec = 0;
minFamSize = 0;
nprod = 0;
%--------------------------------------------------------
srclen = 50; % maximum interaction distance in number of source length between a source and target
dur = 30; % duration of the testing period. 
mthres = 0;
begsttm = datenum(2020,2,1); % starting date of the experiment
endsttm = max(cat0.pricat.datenum); % maximum time reported in the catalog.
sttmlist = [begsttm:dur:endsttm];
FP_ETAS_Times=datetime(sttmlist,'ConvertFrom','datenum');
cd('Outputs/')
save('FP_ETAS_Times.mat','FP_ETAS_Times')
cd ('../')
maxid = length(sttmlist); % id of the testing period. min value for id is 1 and max
%value is N, where N is the maximum number of testing period.
% note that following code can be parallelized if you have access to
% multiple cores.
for i = 1:maxid
    [output] = Wrapper_ETAS_Invert_verCoSeismiq(maxtm,Maux,Mpri,KDE_BKG,...
        OmTyp,btyp,isotyp,prodtyp,srclen,dur,mthres,seqSpec,minFamSize,finiteInteg,i);
end
%% Once all the parameters are inverted use them to simulate catalog

InvUpdateDur = dur; % duration at which inverted parameters are updated.
SimUpdateDur = dur; % durattion at which simulations have to be performed.
nsim = 10^5; % number of simulation.
tresol = dur; % resolution in number of days of the testing period.
maxjobindex = maxid*10; % maximum number of simulation jobs
% note that following code can be parallelized if you have access to
% multiple cores.
for i = 1:maxjobindex
    [output] = Wrapper_Sim_ETAS_verCoSeismiq(nsim,maxtm,InvUpdateDur,SimUpdateDur,...
        Maux,Mpri,KDE_BKG,OmTyp,btyp,isotyp,prodtyp,srclen,mthres,nprod,seqSpec,...
        minFamSz,tresol,i);
end

%% once all the simulations are done perform the evaluation of the models
% we first make a list for the which evaluation has to be performed.
foldnam = './outputvars/';
filepattern2 = 'Eval_sttm%f_Maux%f_Mpri%f_KDEBKG%f_OmTyp%f_btyp%f_isotyp%f_srclen%f_mthres%f_prodtyp%f_seqSpec%f_minFamSz%f_nprod%f_tresol%f.mat';

[missnums] = Wrapper_FailedEval_ver3(foldnam,filepattern2,Maux,Mpri,...
    KDE_BKG,OmTyp,btyp,isotyp,srclen,mthres,prodtyp,nprod,seqSpec,minFamSz,...
    tresol,SimUpdateDur);

save(['EvalList_',Regionnam,'.mat'],'missnums')

%% now perform all the evaluations
LLe_ETAS=zeros(length(FP_ETAS_Times),1); %cumulative log likelihood for each time step
LLe_ETAS_spatial=zeros(1000,length(FP_ETAS_Times)); %1000 cells, hard-coded for now. Log likelihood per cell per time step averages over all simulations
Rate_Events_ETAS=zeros(1000,length(FP_ETAS_Times)); %1000 cells, hard-coded for now. Rate of event per cell per time step averages over all simulations
for i = 1:maxid
    [sim,modeleval] = Wrapper_Eval_ETAS_verCoSeismiq(SimUpdateDur,i);
    LLe_ETAS(i)=modeleval.logliksum;
    LLe_ETAS_spatial(:,i)=modeleval.loglik{1,1}(:,4);
    Rate_Events_ETAS(:,i)=reshape(sim.ratmat{1,1},[1000 1]);
end
cd('Outputs/')
save('LLe_ETAS.mat','LLe_ETAS')
save('LLe_ETAS_spatial.mat','LLe_ETAS_spatial')
save('Rate_Events_ETAS.mat','Rate_Events_ETAS')
cd ('../')




