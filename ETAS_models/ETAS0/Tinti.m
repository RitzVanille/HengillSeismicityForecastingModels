function [beta] = Tinti(mag,M0,dm)
%{
mag: binned magnitude
M0: magnitude of completeness
dm: binsize

beta: bvalue*log(10)
%}
%%
mmag=sum(mag-M0)/length(mag);
beta=1/dm*log(1+dm/mmag);
end