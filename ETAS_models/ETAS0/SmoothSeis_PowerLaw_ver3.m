function [SSM] = SmoothSeis_PowerLaw_ver3(loc,wt,grid,pars,smoothingMethod,lamdbai)
switch smoothingMethod

    case 1
        methodnam = 'Fixed BW';
        SSM = 0; x = round(length(loc(:,1))/100);
        d = 10^pars(1); 
        for i = 1:length(loc(:,1))
            [dist] = HaverSine_fast_ver2(grid,loc(i,:))';
            dist(i) = inf;
            dist = dist.^2;
            kij = (pi^-1) * (pars(2) * d^pars(2)) ./ (dist + d).^(1+pars(2));
            SSM = SSM + wt(i) * kij;
        end
        
       
    case 4
        methodnam = 'Adaptive kernel';
        SSM = 0; x = round(length(loc(:,1))/100);
        d = lamdbai.^2*(10^pars(1));
        
        check = 1;
        
        for i = 1:length(loc(:,1))
            [dist] = HaverSine_fast_ver2(grid,loc(i,:))';
            dist(i) = inf;
            dist = dist.^2;
            kij = (pi^-1) * (pars(2) *d(i)^pars(2)) ./ (dist + d(i)).^(1+pars(2));
            SSM = SSM + wt(i) * kij;
        end
        
end

end