function [curprod,lhat_corr] = UpdateCurProd(lhat,cat,Twin,c,omega,tau)
% first compute the corrected numnber of aftershocks
chekc=1;
F5l      =  Twin(:,1)-(cat(:,1));
F5u      =  Twin(:,2)-(cat(:,1));
F5u(F5u==0) = 1; % minimum time duration
deno      =  (gamma_incomplete_ver2((F5l+c)./tau,-omega)-...
    gamma_incomplete_ver2((F5u+c)./tau,-omega))/...
    gamma_incomplete_ver2(c/tau,-omega);
lhat_corr = lhat./deno;

tempmat = [cat(:,4),lhat_corr];
tempmat = sortrows(tempmat,1);

[umag, iaf] = unique(tempmat(:,1));
[umag, ial] = unique(tempmat(:,1),'last');

curprod = zeros(length(umag),3);
curprod(:,1) = umag;
for i=1:length(umag)
    curprod(i,2) = mean(tempmat(iaf(i):ial(i),2));
    curprod(i,3) = median(tempmat(iaf(i):ial(i),2));
end
chekc  = 1;
lhat_corr = [cat(:,4),lhat_corr];

end