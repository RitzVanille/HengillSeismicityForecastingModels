function [mag]=GR_truncated_simulator_ver2(nsim,Mmax,M0,b)

beta=b*log(10);
n=length(beta);
mag=-log(1+rand(n,nsim).*repmat(-1+exp(-beta.*(Mmax-M0)),[1,nsim]))./...
    repmat(beta,[1,nsim])+repmat(M0,[1,nsim]);

check=1;
end