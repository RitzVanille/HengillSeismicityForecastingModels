function [Tm] = SampleTimeKernel(pars,rnums)


c          = pars(:,1);
om         = pars(:,2);
tau        = 10.^(log10(pars(:,3))-0.1);  % this correction is based on the observation, see Approx_ExpoTaper_PDF.m
A1         = (1./om.*(1./c.^om - 1./(tau+c).^om) + tau./(tau+c).^(1+om)).^-1;
IntegLo    = A1.*1./om.*(1./c.^om - 1./(tau+c).^om);
indlo      = rnums <= IntegLo;
Tm         = zeros(size(rnums));
Tm(indlo)  = (1./c(indlo).^om(indlo) - om(indlo).*rnums(indlo)./A1(indlo)).^(-1./om(indlo)) - c(indlo);
Tm(~indlo) = -tau(~indlo) .* log(1/exp(1) - (rnums(~indlo)-IntegLo(~indlo)).*...
    (tau(~indlo)+c(~indlo)).^(1+om(~indlo))./(A1(~indlo).*tau(~indlo).*exp(1)));


end