function [K] = BranchingRatio2K(n,a,b,Mmax,Mc)
alpha = a/log(10);
tol = 10^-2;

if abs(alpha-b) > tol
    denom = b*(1-10^(-(b-alpha)*(Mmax-Mc)))/((b-alpha)*(1-10^(-b*(Mmax-Mc))));
    
else
    denom = b*log(10)*(Mmax-Mc)/(1-10^(-b*(Mmax-Mc)));
    
end

K = n/denom;
end