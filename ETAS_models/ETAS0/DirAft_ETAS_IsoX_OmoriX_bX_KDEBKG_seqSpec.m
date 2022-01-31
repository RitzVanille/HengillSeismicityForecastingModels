function [DaughtCat] = DirAft_ETAS_IsoX_OmoriX_bX_KDEBKG_seqSpec(MomCat,model,Twin)
%tic

omegatol    =  0.01; % this is to prevent omega from becoming very small and fuck up the simulation
%% assign the ETAS parameters
%if model.seqSpec == 0
bval        =  model.partable.betaGlob/log(10);
M0          =  model.partable.Mpri;
Mmin        =  M0;
Mmax        =  model.Mmax;
%br          =  model.partable.n;
K           =  10.^model.partable.KNorm;
a           =  model.partable.aNorm;
d           =  10.^model.partable.d;
gam         =  model.partable.gamma;
%if model.isotyp == 1 || model.isotyp == 2 || model.isotyp == 3
rho         =  model.partable.rho;
dm = model.partable.dm;


magbinid = round((MomCat(:,4) - model.Maux)/model.dm)+1;
umags = BinMags([model.Maux:model.dm:11]',0,model.dm);

switch model.OmTyp
    
    case 1
        tau   =  repmat(10.^model.partable.tau0,size(MomCat(:,4)));
        c     =  repmat(10.^model.partable.c0,size(MomCat(:,4)));
        omega =  repmat(model.partable.p0-1,size(MomCat(:,4)));
   
end

omega(abs(omega)<0.01) = sign(omega(abs(omega)<0.01))*omegatol;

%% find the expected number of triggered events
% if model.Maux~=model.M0

%F9   =  MomCat(:,4)-Mmin;
F9   =  umags-Mmin;


    totaftnum = K.*exp(a.*F9);
totaftnum1 = totaftnum(magbinid);

numaft = poissrnd(totaftnum1);
%toc
numaft(:,2)              = [1:length(numaft(:,1))]';
numaft(numaft(:,1)==0,:) = [];
n                        = sum(numaft(:,1));
%%
%toc
if n ~= 0
    repmainid1 = cumsum(numaft(:,1));
    repmainid2 = [0;repmainid1(1:end-1)]+1;
    iafl       = [repmainid2,repmainid1];
    mainid     = zeros(n,1);
    for i = 1:size(iafl,1)
        mainid(iafl(i,1):iafl(i,2),1) = numaft(i,2);
    end
    
    c     = c(mainid);
    omega = omega(mainid);
    tau   = tau(mainid);
    ompars = [c,omega,tau];
    
    %% simulate time
    randnum_om = rand(n,1);
    [afttm]    = SampleTimeKernel_ver5(ompars,randnum_om,model.OmTyp);
    afttm(:,2) = mainid;
    afttm(:,1) = afttm(:,1)+MomCat(afttm(:,2),1);
    rmind      = afttm(:,1)<Twin(1) | afttm(:,1)>Twin(2);
    afttm(rmind,:) = [];
    
    
    
    if ~isempty(afttm)
        
        n = size(afttm,1);
        
        
        %if (model.isotyp == 1 || model.isotyp == 2 || model.isotyp == 4) && model.simtyp == 1 && ...
        %        (model.prodtyp == 0 || model.prodtyp == 1 || model.prodtyp == 3 || model.prodtyp == 301 || model.prodtyp == 302)
        
        if model.isotyp == 4
            rho = rho(afttm(:,2));
        end
        %% simulate location
        randnum_loc   = rand(size(afttm,1),1);
        randnum_theta = rand(size(afttm,1),1)*360;
        r             = sqrt(d*exp(gam*(MomCat(afttm(:,2),4)-M0)).*(1./(1-randnum_loc).^(1./rho)-1));
        aftcat        = [afttm(:,1),r,MomCat(afttm(:,2),[2,3])];
        degdist       = aftcat(:,2)*180/(6371*pi);
        [latout,lonout] = reckon(aftcat(:,4),aftcat(:,3),degdist,randnum_theta);
        
        %% simulate magnitude
        chec = 1;
        
        switch model.btyp
            case 1
                nbval = repmat(bval(1),[n,1]);
                [mag] = GR_truncated_simulator_ver2(1,repmat(Mmax,[n,1]),...
                    repmat(M0-dm/2,[n,1]),nbval);
  
        end
        
        [mag] = BinMags(mag,0,dm);
        
        %%
        DaughtCat = [aftcat(:,1),lonout,latout,mag,zeros(size(mag)),zeros(size(mag)),ones(size(mag)),MomCat(afttm(:,2),8)];
    else
        DaughtCat = [];
    end
else
    DaughtCat = [];
end

end