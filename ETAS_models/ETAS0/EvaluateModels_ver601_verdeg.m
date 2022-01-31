function [loglik,logliksum] = EvaluateModels_ver601_verdeg(sim,div)
% smoothing is done using power law kernel
%%
logliksum = zeros(length(sim.pixlev),length(sim.Mt));
loglik    = cell(length(sim.pixlev),length(sim.Mt));

%tic
for i = 1:length(sim.pixlev)
    %i
    for j = 1:length(sim.Mt)
        
        disp(strcat(['PixLev : ', num2str(sim.pixlev(i)), ...
            '; Mt : ',num2str(sim.Mt(j))]))
        
        
        AllId      =    sim.tformgrid{i};
        urealtrid  =    sim.ValidCount{i,j}(:,1:2);
        [Lia,Locb] =    ismember(AllId,urealtrid,'rows');
        obscount   =    sim.ValidCount{i,j}(:,3);
        maxpossnum   =    sum(obscount)+1;
        
        AllId(:,3)   =    0;
        AllId(Lia,3) =    obscount;
        templik      =    zeros(length(AllId(:,1)),1);
        
        [~,id] =    ismember(sim.SimCount{i,j}(:,1:2),AllId(:,1:2),'rows');
        
        for k = 1:length(AllId(:,1))
            
            tselempprob = sim.SimCount{i,j}(id == k,:);
            
            
            if ~isempty(tselempprob)
                
                selempprob = [0,sim.nsim-sum(tselempprob(:,3))+1;tselempprob(:,3:4)];
                
            else
                selempprob = [0,sim.nsim-1/1000;1,1/1000];
            end
            %selempprob(selempprob(:,3)==0,:) = [];

            nsupport   =  [0:1:max(maxpossnum,max(selempprob(:,1)))]';
            obsind = selempprob(:,1) == AllId(k,3);
            % this if statement makes things faster, as we dont have to do
            % any smoothing
            if sum(obsind) == 1
                
                templik(k) = log(selempprob(obsind,2)/sim.nsim);
            else
                emprob            = zeros(size(nsupport));
                indtyp1           = selempprob(selempprob(:,2)>0,1);
                emprob(indtyp1+1) = selempprob(selempprob(:,2)>0,2)/sim.nsim;
                indtyp2           = selempprob(selempprob(:,2)<=0,1);
                tvec              = [min(selempprob(:,1)):1:max(selempprob(:,1))]';
                holid             = tvec(~ismember(tvec,selempprob(:,1)));
                indtyp2           = [indtyp2;holid];
                
                if ~isempty(indtyp2)
                    support   =  [-0.5,max(maxpossnum,max(selempprob(:,1)))];
                    sig       =  ones(size(selempprob,1),1);
                    for k1 = 1:length(selempprob(:,1))
                        sig(k1) = sum(abs(selempprob(:,1)-selempprob(k1,1)).*...
                            selempprob(:,2))/sim.nsim;
                    end
                    sig  =  sig/div;
                    
                    
                    allnum    =    indtyp2;
                    smpdf     =    zeros(length(indtyp2),1);
                    for k2  =  1:length(allnum)
                        ny      =  0;
                        numobs  =  allnum(k2);
                        for k3  =  1:size(selempprob,1)
                            
                            y   =  normcdf(numobs+0.5,selempprob(k3,1),sig(k3))-...
                                normcdf(numobs-0.5,selempprob(k3,1),sig(k3));
                            y   =  selempprob(k3,2).*y/...
                                (normcdf(support(end),selempprob(k3,1),sig(k3))-...
                                normcdf(support(1),selempprob(k3,1),sig(k3)));
                            ny  =   ny+y;
                            
                        end
                        ny        =  ny/sim.nsim;
                        smpdf(k2)  =  ny;
                    end
                    emprob(indtyp2+1)  =  smpdf;
                end
                emprob  =    emprob/sum(emprob);
                logempProb = zeros(size(emprob));
                
                %meanrat  =   sum(selempprob(:,2).*selempprob(:,3))./sim.nsim;
                indtyp3  =   find(emprob==0);
                Notindtyp3  =   find(emprob~=0);
                logempProb  = log(emprob);
                Idx         =   knnsearch(Notindtyp3,indtyp3,'K',1);
                % extrapolate using poissonian from the closet bin
                if ~isempty(indtyp3)
                    poissExtrap = zeros(size(indtyp3));
                    for k4 = 1:length(indtyp3)
                        x = (indtyp3(k4)-1);
                        lambda = Idx(k4)-1;
                        
                        poissExtrap(k4) = (x*log(lambda)-lambda - ...
                            sum(log(1:x))+log(emprob(Idx(k4))));
                    end
                    logempProb(indtyp3) = poissExtrap;
                end
                nnorm = log(sum(exp(logempProb)));
                logempProb = logempProb-nnorm;
                check = 1;
                templik(k) = logempProb(AllId(k,3)+1);
                
            end
        end
        totrumtm        = toc;
        disp(strcat('Runtime : ' , num2str(totrumtm)))
        AllId(:,4)      = templik;
        logliksum(i,j)  = sum(templik);
        disp(strcat('LL : ' , num2str(sum(templik))))
        disp('----------------------------')
        loglik{i,j}     = AllId;
    end
    
    
end


%%

end
