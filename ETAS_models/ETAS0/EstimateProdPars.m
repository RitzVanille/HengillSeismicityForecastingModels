function [log10K,a] = EstimateProdPars(lhat,M0)

[~,Krang,arang] = ETAS_Globe_Par_Range('narrow');
initpars        = [rand*(Krang(2)-Krang(1))+Krang(1),...
    rand*(arang(2)-arang(1))+arang(1)];

[~,Krang,arang] = ETAS_Globe_Par_Range('wide');
nrang           = [Krang',arang'];
FixedQuant2     = -(gammaln(lhat(:,2)+1));
options = optimoptions('fmincon','Algorithm','interior-point','Display','iter',...
    'TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',100*length(initpars),'MaxIter',100);
[fpars,NLL] = ...
    fmincon(@ProdLL,initpars,[],[],[],[],nrang(1,:),nrang(2,:),[],options);

log10K = fpars(1);
a      = fpars(2);

%%
    function [NLL] = ProdLL(pars)
        K0_f1     = 10.^pars(1);
        a_f1      = pars(2);
        G         = K0_f1.*exp(a_f1.*(lhat(:,1)-M0));
        NLL       = -sum(FixedQuant2 -G + lhat(:,2).*log(G));
    end

end