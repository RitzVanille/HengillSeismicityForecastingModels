function [SRL,SSRL,RW,RA,AD] = coppersmith_ver2(mag,type)

switch type
    case 1
        % this the case of strike slip faults
        % surface rupture length
        SRL=10.^(0.74*mag-3.55);
        % subsurface rupture length
        SSRL=10.^(0.62*mag-2.57);
        % rupture width
        RW=10.^(0.27*mag-0.76);
        % rupture area
        RA=10.^(0.9*mag-3.42);
        % average slip
        AD=10.^(0.9*mag-6.32);
    case 2
        % this is for reverse faults
        % surface rupture length
        SRL=10.^(0.63*mag-2.86);
        % subsurface rupture length
        SSRL=10.^(0.58*mag-2.42);
        % rupture width
        RW=10.^(0.41*mag-1.61);
        % rupture area
        RA=10.^(0.98*mag-3.99);
        % average slip
        AD=10.^(0.08*mag-0.74);
    case 3
        % this is for normal faults
        % surface rupture length
        SRL=10.^(0.5*mag-2.01);
        % subsurface rupture length
        SSRL=10.^(0.5*mag-1.88);
        % rupture width
        RW=10.^(0.35*mag-1.14);
        % rupture area
        RA=10.^(0.82*mag-2.87);
        % average slip
        AD=10.^(0.63*mag-4.45);
    case 4
        % this is for oblique faults
        % surface rupture length
        SRL=10.^(0.69*mag-3.22);
        % subsurface rupture length
        SSRL=10.^(0.59*mag-2.44);
        % rupture width
        RW=10.^(0.32*mag-1.01);
        % rupture area
        RA=10.^(0.91*mag-3.49);
        % average slip
        AD=10.^(0.69*mag-4.80);
end

    
end