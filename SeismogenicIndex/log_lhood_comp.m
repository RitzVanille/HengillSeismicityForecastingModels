function log_lhood_comp = log_lhood_comp(theta,rate)

af = theta(1);
%tau = theta(2);
b = theta(2); % if relaxation time tau acivated: b=theta(3)

%
A1 = rate.N*(af-b*rate.m_0)/log10(exp(1));          % first term in the likelihood 
%
dotV_bs_ts = interp1(rate.t_b_s,rate.dot_V_bs,rate.t_sbs);
K2 = sum(log(dotV_bs_ts));              % second term in the likelihood
if K2<1.e30
    dotV_bs_ts=dotV_bs_ts+1e-10;
    K2 = sum(log(dotV_bs_ts));
end
    
%
K3 = 0; %(rate.N-rate.N_sbs)*log(rate.dot_V_shut_in);
%
A4 = 0; %-1/tau*(sum(rate.t_sas-rate.T_s));
%
A5 = -10^(af-b*rate.m_0).*rate.tot_V; %+rate.dot_V_shut_in*tau*(1-exp(-(rate.T-rate.T_s)/tau)));
%
A6 = rate.N*log(b);
%
K4 = rate.N*log(log(10));
%
A7 = -b*log(10)*sum(rate.data_magn);
%
A8 = rate.N*b*rate.m_0_m*log(10);
% Magnitude terms 
log_lhood_comp =  -(A1+K2+K3+A4+A5+A6+K4+A7+A8);

end