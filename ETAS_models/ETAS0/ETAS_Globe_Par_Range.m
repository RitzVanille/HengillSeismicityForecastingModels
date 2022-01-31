function [murang,Krang,arang,drang,rorang,gammarang,...
    c0rang,c1rang,p0rang,p1rang,taurang,betabkgrang,...
    betalrang,betaurang,tau0rang,tau1rang,alpha_slope_rang,...
    alpha_inter_rang,q0rang,q1rang,M0rang,Delrang,p2rang...
    ,alpha0rang,alpha1rang,alpha2rang] = ETAS_Globe_Par_Range(rangetyp)

switch rangetyp
    case 'narrow'
        murang     = [-12,-8];
        Krang      = [-1,0];
        arang      = [2,3];
        
        drang      = [1,4];
        rorang     = [0.5,2];
        gammarang  = [1,2];
        
        c0rang     = [-4,-2];
        c1rang     = [-2,3];
        p0rang     = [1.01,1.5];
        p1rang     = [-1,1];
        taurang    = [3,4];
        tau0rang    = [3,4];
        tau1rang    = [0,1];
        
        betabkgrang = [2, 3];
        betalrang   = [1.5, 2];
        betaurang   = [2.5, 3];
        alpha_slope_rang     = [0.2,0.6];
        alpha_inter_rang     = [-2,-1];
        q0rang     = [1.01,2.5];
        q1rang     = [-1,1];
        
        M0rang = [1,2.5];
        Delrang = [0,1];
        
        p2rang = [-1,1];
        
        alpha0rang = [1,3];
        alpha1rang = [0,1];
        alpha2rang = [0,1];
    case 'wide'
        murang      = [-12,3];
        Krang       = [-12,5];
        arang       = [0,10];
        
        drang       = [-3,6];
        rorang      = [0,5];
        gammarang   = [0,5];
        
        c0rang      = [-13,1];
        c1rang      = [-3,3];
        p0rang      = [-1,2];
        p1rang      = [-2,2];
        taurang     = [0,6];
        tau0rang    = [0,6];
        tau1rang    = [-1,3];
        
        betabkgrang = [0.1, 5];
        betalrang   = [0.1, 5];
        betaurang   = [0.1, 5];
        alpha_slope_rang     = [0,5];
        alpha_inter_rang     = [-10,10];
        q0rang      = [-6,6];
        q1rang      = [-2,4];
        
        M0rang = [0,4];
        Delrang = [0,2];
        p2rang = [-2,2];
        
        alpha0rang = [0,20];
        alpha1rang = [0,2];
        alpha2rang = [0,2];
    case 'timeonly_narrow'
        murang     = [-3,0];
        Krang      = [-1,0];
        arang      = [2,3];
        
        drang      = [1,3];
        rorang     = [0.5,2];
        gammarang  = [1,2];
        
        c0rang     = [-4,-2];
        c1rang     = [-2,3];
        p0rang     = [1.01,1.5];
        p1rang     = [-1,1];
        taurang    = [4,6];
        tau0rang    = [3,4];
        tau1rang    = [0,1];
        
        betabkgrang = [2, 3];
        betalrang   = [1.5, 2];
        betaurang   = [2.5, 3];
end


end