function [fpars,LL,methodnam,brkflag] = InvertSmoothingPar_PowerLaw_ver3(wtsrc,wtrec,globvars,method,hyperpars,initpars,maxrad,lambdai)
globvars(globvars(:,1) == 0, : ) = [];
globwt = wtsrc(globvars(:,2));
globvars(globwt==0,:) = [];
globwt(globwt==0) = [];
[Cglob,iaf]  = unique(globvars(:,3));
[Cglob,ial]  = unique(globvars(:,3),'last');


globmaxrad = maxrad(globvars(:,2));
%globmaxrad =
recwt = wtrec(Cglob);
Norm = sum(wtsrc);
switch method
    case 1
        
        brkflag = 0;
        methodnam = 'Fixed BW';
        mind = hyperpars(1);
        drang = log10([mind^2,5000]);
        rorang = [0, 2];
        if isempty(initpars)
            initpars = [rand*(drang(2)-drang(1))+drang(1),...
                rand*(rorang(2)-rorang(1))+rorang(1)];
        end
        drang = log10([mind^2,inf]);
        rorang = [0, 10];
        nrang = [drang',rorang'];
        options = optimoptions('fmincon','Algorithm','interior-point','Display','off',...
            'TolFun',1e-6,'TolX',1e-6,'MaxFunEvals',100*length(initpars),'MaxIter',100);
        [fpars,LL] = ...
            fmincon(@OptimFixedBW,initpars,[],[],[],[],nrang(1,:),nrang(2,:),[],options);
        check = 1;
        
        
    case 4
        
        brkflag = 0;
        methodnam = 'Adaptive kernel';
        %
        globlambda = lambdai(globvars(:,2)).^2;
        mind = hyperpars(1);
        drang = log10([mind^2,5000]);
        rorang = [0, 2];
        if isempty(initpars)
            initpars = [rand*(drang(2)-drang(1))+drang(1),...
                rand*(rorang(2)-rorang(1))+rorang(1)];
        end
        drang = log10([mind^2,inf]);
        rorang = [0, 10];
        nrang = [drang',rorang'];
        options = optimoptions('fmincon','Algorithm','interior-point','Display','iter',...
            'TolFun',1e-6,'TolX',1e-6,'MaxFunEvals',100*length(initpars),'MaxIter',100);
        initpars
        [fpars,LL] = ...
            fmincon(@OptimAdaptKernel,initpars,[],[],[],[],nrang(1,:),nrang(2,:),[],options);
        chekc = 1;
        
end

check=1;

%%
    function [LL] = OptimFixedBW(pars)
        d = 10^pars(1);
        ro = pars(2);
        %normnew = (pi^-1) * (ro * d^ro);
        normnew = ro/pi*1./(1./ d.^ro - 1./(globmaxrad+d).^ro);
        kij = normnew ./ (globvars(:,1) + d).^(1+ro);
        PDFij = globwt .* kij;
        SPDF = nan(size(Cglob));
        for i = 1:length(Cglob)
            SPDF(i) = sum(PDFij(iaf(i):ial(i)));
        end
        LL = sum(recwt .*  log(SPDF)) - sum(recwt) * log(Norm);
        LL = -LL;
        check=1;
    end
%%
    function [LL] = OptimAdaptKernel(pars)
        h = 10^pars(1);
        ro = pars(2);
        
        d = h * globlambda;
        %d(d<100) = 100;
        %normnew = (pi^-1) * (ro * d.^ro);
        normnew = ro/pi*1./(1./ d.^ro - 1./(globmaxrad+d).^ro);
        kij = normnew ./ (globvars(:,1) + d).^(1+ro);
        PDFij = globwt .* kij;
        SPDF = nan(size(Cglob));
        for i = 1:length(Cglob)
            SPDF(i) = sum(PDFij(iaf(i):ial(i)));
        end
        LL = sum(recwt .*  log(SPDF)) - sum(recwt) * log(Norm);
        LL = -LL;
        check=1;
        
    end

end