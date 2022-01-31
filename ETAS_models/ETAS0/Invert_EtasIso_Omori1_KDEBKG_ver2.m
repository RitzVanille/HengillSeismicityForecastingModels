function [pars,LLall,trigprob,lhat_estep,IP,prodTable,lhat_corr,ProdPars,SMpars,Adaptive_BW,mu_svar,Invcomp] =...
    Invert_EtasIso_Omori1_KDEBKG_ver2(cat,M0,globvars,Avor1,timewin,catrec,initpars,initmodel,smHyper,globoptions,globvars_bkg,maxrad)
M1=M0;
dm=0.1;
rho_prior = 0.5;
alpha_sens = 0.5;
mind = smHyper(3);
d_prior = 10;
minAftBuffer = 7;

%% obtain the initial ssm
wtsrc = ones(length(catrec(:,1)),1)*0.5;
wtrec = wtsrc;
initSMpars = [log10(d_prior^2),rho_prior];
if ~isempty(initmodel)
    if ~isempty(initmodel.IP)
        wtsrc = ones(length(catrec(:,1)),1)*0.5;
        aptid = min([length(wtsrc),length(initmodel.IP)]);
        wtsrc(1:aptid) = initmodel.IP(1:aptid);
        wtrec = wtsrc;
    end
    if ~isempty(initmodel.SMpars)
        initSMpars = initmodel.SMpars;
    end
end



[SMpars,LL_SSM,~,brkflag] = InvertSmoothingPar_PowerLaw_ver3(wtsrc,wtrec,...
    globvars_bkg,1,mind,initSMpars,maxrad,[]);
initSMpars = SMpars;
old_LL_SSM = LL_SSM;
[mu_svar] = SmoothSeis_PowerLaw_ver3(catrec(:,2:3),wtsrc,catrec(:,2:3),SMpars,1,[]);
mu_svar = mu_svar/diff(timewin);

%% define the initial value and range of parameters
% first define intial value and range of all parameters except mu
indcatrec =  (ismember(cat,catrec,'rows'));
ntwin=[max(cat(:,1),timewin(1,1)),...
    repmat(timewin(1,2),length(cat(:,1)),1)]; % cat is the source catalog
indbufferok = ntwin(:,2) - ntwin(:,1) >= minAftBuffer;
ntwin = ntwin(indcatrec & indbufferok,:);
selcatsrc = cat(indcatrec & indbufferok,:);
reliableid = 1;


tmlen=timewin(2)-timewin(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[murang,Krang,arang,drang,rorang,gammarang,...
    crang,~,p0rang,~,taurang] = ETAS_Globe_Par_Range('narrow');

if isempty(initpars)
    initpars=[rand*(murang(2)-murang(1))+murang(1),...
        rand*(Krang(2)-Krang(1))+Krang(1),...
        rand*(arang(2)-arang(1))+arang(1),...
        rand*(drang(2)-drang(1))+drang(1),...
        rand*(rorang(2)-rorang(1))+rorang(1),...
        rand*(gammarang(2)-gammarang(1))+gammarang(1),...
        rand*(crang(2)-crang(1))+crang(1),...
        rand*(p0rang(2)-p0rang(1))+p0rang(1),...
        rand*(taurang(2)-taurang(1))+taurang(1)];
end


[~,Krang,arang,drang,rorang,gammarang,...
    crang,~,p0rang,~,taurang] = ETAS_Globe_Par_Range('wide');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Cglob,iafglob]=unique(globvars(:,4));
[Cglob,ialglob]=unique(globvars(:,4),'last');



%% starting inversion
parsold=initpars;
parsnew=inf;
tolpars=10^-2;
tolLL=1;
diffLL_f=inf;
difpars=inf;
kfiter=1;
%trigthres=10^-50;
disp('==============Begining EM ========================')
[trigprob,lhat_estep,nhat_estep,IP]=Estep(parsold);

disp('Smoothing new IPs')
wtsrc = IP;
wtrec = wtsrc;

[SMpars,LL_SSM,~,brkflag] = InvertSmoothingPar_PowerLaw_ver3(wtsrc,wtrec,...
    globvars_bkg,1,mind,initSMpars,maxrad,[]);
initSMpars = SMpars;
old_LL_SSM = LL_SSM;
[mu_svar] = SmoothSeis_PowerLaw_ver3(catrec(:,2:3),wtsrc,catrec(:,2:3),SMpars,1,[]);
mu_svar = mu_svar/diff(timewin);

nhat_estep(nhat_estep<=0)=10^-85;
mupars=log10(nhat_estep./(Avor1*tmlen));
trigprob=[trigprob,globvars];
%trigprob(trigprob(:,1)<trigthres,:)=[];
disp('===================================================')
old_LL_f=inf;
Invcomp = 1;
while (difpars>globoptions.tolpars || diffLL_f>globoptions.tolLL) && kfiter <= globoptions.maxiter ...
        && toc(globoptions.tmpars(1))<= globoptions.tmpars(2)
    %disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    %tic
    [ETAS_pars,LL_f]=Mstep(parsold(2:end));
    curLL = LL_f + LL_SSM;
    parsoldcopy=parsold;
    parsold = [mupars,ETAS_pars];
    
    parsnew = parsold;
    difpars = sum(abs(parsoldcopy-parsnew));
    
    disp(strcat('Following are results from EM step::::',num2str(kfiter)))
    
    disp('=======================')
    disp('Accuracy::::')
    accuracy = log10(difpars);
    disp(accuracy)
    parsold = parsnew;
    pars4display = parsnew;
    disp('=======================')
    disp('ETAS pars::::')
    disp(pars4display)
    %LL_f=LL_final+LL_final1+LL_final2+LL_final3+LL_final4;
    diffLL_f = log10(abs(curLL-old_LL_f));
    disp('=======================')
    disp('LL::::')
    disp(diffLL_f)
    old_LL_f = curLL
    
    %TimeTaken(kfiter)=toc;
    kfiter=kfiter+1;
    [trigprob,lhat_estep,nhat_estep,IP]=Estep(parsold);
    disp('Smoothing new IPs')
    wtsrc = IP;
    wtrec = wtsrc;
    
    [SMpars,LL_SSM,~,brkflag] = InvertSmoothingPar_PowerLaw_ver3(wtsrc,wtrec,...
        globvars_bkg,1,mind,initSMpars,maxrad,[]);
    initSMpars = SMpars;
    old_LL_SSM = LL_SSM;
    [mu_svar] = SmoothSeis_PowerLaw_ver3(catrec(:,2:3),wtsrc,catrec(:,2:3),SMpars,1,[]);
    mu_svar = mu_svar/diff(timewin);
    disp('=======================')
    SMpars
    disp('=======================')
    nhat_estep(nhat_estep<=0)=10^-85;
    disp(strcat('current branching ratio : ',num2str(1-sum(IP)/length(IP))));
    disp(strcat('current background number : ',num2str(sum(IP))));
    disp('===================================================================')
    mupars=log10(nhat_estep./(Avor1*tmlen));
    trigprob=[trigprob,globvars];
    if toc(globoptions.tmpars(1)) > globoptions.tmpars(2)
        Invcomp = 0;
    end
    
    %     if difpars<globoptions.tolpars
    %         disp('Nearly no change in parameters since the last step')
    %         break
    %     end
    %trigprob(trigprob(:,1)<trigthres,:)=[];
    %toc
end
disp('convergence criteria satified...')

disp('======================END EM ========================')
pars=[parsnew];

Adaptive_BW = [];

[prodTable,lhat_corr] = UpdateCurProd(lhat_estep,selcatsrc,ntwin,...
    10^pars(7),pars(8)-1,10^pars(9));
[log10K,a] = EstimateProdPars(lhat_corr,M0);
ProdPars   = [log10K,a];
disp(ProdPars)

tp1 = zeros(size(nhat_estep));
tp1(nhat_estep>0) = nhat_estep(nhat_estep>0).*log(nhat_estep(nhat_estep>0));

LLbkg = -sum(-(gammaln(nhat_estep+1))-nhat_estep+tp1);

LLall = LL_f+LLbkg;
LLall= [-LLall, -LL_SSM];

check=1;
disp('Returning the output....')
disp(' :) Hopefully they are good :) ')
disp('===================================================')




%% Estep
    function [trigprob_Estep,lhat_Estep,nhat_Estep,IP]=Estep(parsEstep)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dEstep=10^parsEstep(4);  % in log
        roEstep=parsEstep(5);
        gammaEstep=parsEstep(6);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        muEstep = mu_svar ;
        K0Estep=10.^parsEstep(2);
        aEstep=parsEstep(3);
        cEstep=10.^parsEstep(7);
        omEstep=parsEstep(8)-1;
        tauEstep=10.^parsEstep(9);
        
        
        numer=K0Estep.*(exp(-globvars(:,2)./tauEstep)./...
            (globvars(:,2)+cEstep).^(1+omEstep)).*...
            (1./(globvars(:,1)+dEstep*...
            exp(gammaEstep*(globvars(:,3)-M0))).^(1+roEstep)).*...
            exp(aEstep.*(globvars(:,3)-M0));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % computing triggering probability
        kin=1;
        kfin=1;
        deno=zeros(size(numer,1),1);
        for i=1:length(Cglob)
            deno(iafglob(i):ialglob(i),1)=sum(numer(iafglob(i):ialglob(i)))...
                +muEstep(Cglob(i));
            
        end
        trigprob_Estep=numer./deno(:,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%% computing lhat
        %lhat_Estep is the number of aftershocks trigerred by earthquakes
        %it the source catalog
        trigprob_Estep1=[trigprob_Estep,globvars(:,5)];
        %trigprob_Estep1(trigprob_Estep<trigthres,:)=[];
        trigprob_Estep1=sortrows(trigprob_Estep1,2);
        [~,iaf]=unique(trigprob_Estep1(:,2));
        [Cind,ial]=unique(trigprob_Estep1(:,2),'last');
        lhat_Estep=zeros(length(cat(:,1)),1);
        for s1=1:length(Cind)
            lhat_Estep(Cind(s1),1)=sum(trigprob_Estep1(iaf(s1):ial(s1),1));
        end
        
        lhat_Estep = lhat_Estep(indcatrec & indbufferok);
        %nhat_Estep=-sum(lhat_Estep)+size(cat,1);
        
        %%%%%%%%%%% computing independence probability
        trigprob_Estep1=[trigprob_Estep,globvars(:,4)];
        [~,iaf]=unique(trigprob_Estep1(:,2));
        [Cind,ial]=unique(trigprob_Estep1(:,2),'last');
        
        IP=nan(size(catrec,1),1);
        for i2=1:length(Cind)
            IP(Cind(i2))=1-sum(trigprob_Estep1(iaf(i2):ial(i2),1));
        end
        IP(isnan(IP))=1;
        IP(IP<0)=0;
        nhat_Estep=sum(IP);
        check1=1;
    end


%%
    function [finalpars_Mstep,LL_Mstep]=Mstep(oldpars_Mstep)
        nrang=[Krang',arang',drang',rorang',gammarang',crang',p0rang',taurang'];
        ninit=oldpars_Mstep;
        options = globoptions.options;
        % Option to display output
        FixedSum2=sum(trigprob(:,1).*(trigprob(:,4)-M0));
        FixedSum1=sum(trigprob(:,1));
        FixedQuant1=(trigprob(:,4)-M0);
        
        lhat_f=lhat_estep;
        
        PC2=sum(lhat_f);
        PC1=sum(trigprob(:,1).*trigprob(:,3));
        
        F5l=ntwin(:,1)-(selcatsrc(:,1));
        F5u=ntwin(:,2)-(selcatsrc(:,1));
        F8=selcatsrc(:,4)-M0;
        FixedQuant2=-(gammaln(lhat_f+1));
        
        [finalpars_Mstep,LL_Mstep] = ...
            fmincon(@loglik_Mstep,ninit,[],[],[],[],nrang(1,:),nrang(2,:),[],options);
        function [LL]=loglik_Mstep(pars_f)
            
            
            K0_f1=10.^pars_f(1);
            a_f1=pars_f(2);
            d_f=10^pars_f(3);  % in log
            ro_f=pars_f(4);
            gamma_f=pars_f(5);
            c_f=10^pars_f(6);  % in log
            om_f=pars_f(7)-1;
            tau_f=10^pars_f(8);
            
            
            
            G=K0_f1.*pi.*exp(a_f1.*F8)./...
                (ro_f*((d_f*exp(gamma_f*F8)).^ro_f));
            
            %G=K0_f1.*exp(a_f1.*F8).*F9;
            
            G=G.*exp(c_f./tau_f).*(tau_f.^(-om_f)).*...
                (gamma_incomplete_ver2((F5l+c_f)./tau_f,-om_f)-...
                gamma_incomplete_ver2((F5u+c_f)./tau_f,-om_f));
            G(G<=0)=10^-100;
            term2=sum(FixedQuant2-G+lhat_f.*log(G));
            
            term3=FixedSum1*(log(ro_f)+ro_f*log(d_f)-log(pi))-...
                sum(trigprob(:,1).*log((trigprob(:,2)+d_f*exp(gamma_f*FixedQuant1))))*(1+ro_f)+...
                FixedSum2*ro_f*gamma_f;
            
            term1=-sum(trigprob(:,1).*log((trigprob(:,3)+c_f)))*(1+om_f)-PC1/tau_f-...
                (c_f/tau_f-(om_f)*log(tau_f)+log(gamma_incomplete_ver2(c_f/tau_f,-om_f)))*PC2;
            
            LL=term1+term2+term3;
            LL=-LL;
        end
        
    end
end
