function [p, y0, tspan, options] = MultipleFcRParameters(type1,type2,disease,bodycell)
k1af = 10e-6; %IgG1 + Ag forward
k1ar = 2e-4;
k2af = 10e-6; %IgG2 + Ag forward
k2ar = 2e-4;
k3af = 10e-6; %IgG3 + Ag forward
k3ar = 2e-4;
k4af = 10e-6; %IgG4 + Ag forward
k4ar = 2e-4;

% k1r2f = 52e-6; %IgG1:Ag + FcR2a forward
k1r2r = 1e-2;
% k2r2f = 4.5e-6; %IgG2:Ag + FcR2a forward
k2r2r = 1e-2;
% k3r2f = 8.9e-6; %IgG3:Ag + FcR2a forward
k3r2r = 1e-2;
% k4r2f = 1.7e-6; %IgG1:Ag + FcR2a forward
k4r2r = 1e-2;

% k1r3f = 20e-6; %IgG1:Ag + FcR3a forward
k1r3r = 1e-2;
% k2r3f = .7e-6; %IgG2:Ag + FcR3a forward
k2r3r = 1e-2;
% k3r3f = 98e-6; %IgG3:Ag + FcR3a forward
k3r3r = 1e-2;
% k4r3f = 2.5e-6; %IgG1:Ag + FcR3a forward
k4r3r = 1e-2;

if type1 == '2aH'
    k1r2f = 52e-6;
    k2r2f = 4.5e-6;
    k3r2f = 8.9e-6;
    k4r2f = 1.7e-6;    
elseif type1 == '2aR'
    k1r2f = 35e-6;
    k2r2f = 1e-6;
    k3r2f = 9.1e-6;
    k4r2f = 2.1e-6;
elseif type1 == '3aV'
    k1r2f = 20e-6;
    k2r2f = .7e-6;
    k3r2f = 98e-6;
    k4r2f = 2.5e-6;
elseif type1 == '3aF'
    k1r2f = 11.7e-6;
    k2r2f = .3e-6;
    k3r2f = 77e-6;
    k4r2f = 2e-6;
elseif type1 == 'fcr1'
    k1r2f = 6e-4;
    k2r2f = 0;
    k3r2f = 6e-4;
    k4r2f = 3e-4;

    k2r2r = 0;
end

if type2 == '2aH'
    k1r3f = 52e-6;
    k2r3f = 4.5e-6;
    k3r3f = 8.9e-6;
    k4r3f = 1.7e-6;    
elseif type2 == '2aR'
    k1r3f = 35e-6;
    k2r3f = 1e-6;
    k3r3f = 9.1e-6;
    k4r3f = 2.1e-6;
elseif type2 == '3aV'
    k1r3f = 20e-6;
    k2r3f = .7e-6;
    k3r3f = 98e-6;
    k4r3f = 2.5e-6;
elseif type2 == '3aF'
    k1r3f = 11.7e-6;
    k2r3f = .3e-6;
    k3r3f = 77e-6;
    k4r3f = 2e-6;
elseif type2 == 'fcr1'
    k1r3f = 6e-4;
    k2r3f = 0;
    k3r3f = 6e-4;
    k4r3f = 3e-4;

    k2r3r = 0;
end

% new total concentrations for IgGs from average of converted MFI values
if disease == 'Sjogrens'
    if bodycell == 'assay'
        igg1tot = 9.46;
        igg2tot = 1.8693;
        igg3tot = 1.4294;
        igg4tot = .2387;
        agtot = 12;
        fcr2atot = 8.67/2;
        fcr3atot = 8.67/2;
    elseif bodycell == 'blood'
        igg1tot = 9.46;
        igg2tot = 1.8693;
        igg3tot = 1.4294;
        igg4tot = .2387;
        agtot = 12;
        fcr2atot = 8.67/2;
        fcr3atot = 8.67/2;
    end
elseif disease == 'HIV'
    % total concentrations for HIV
    if bodycell == 'assay'
        igg1tot = 85.4195;
        igg2tot = 0.132837;
        igg3tot = 3.06998;
        igg4tot = 0.102194;
        agtot = 25;
        fcr2atot = 20/2;
        fcr3atot = 20/2;
    elseif bodycell == 'blood'
        igg1tot = 581;
        igg2tot = 2.3;
        igg3tot = 82.6;
        igg4tot = 0.23;
        agtot = 2.2e-6;
        fcr2atot = 2.5e-2;
        fcr3atot = 5.5e-2;
    elseif bodycell == 'rectal'
        igg1tot = 0.277;
        igg2tot = 0.0011;
        igg3tot = 0.092;
        igg4tot = 0.00011;
        agtot = 1.1e-5;
        fcr2atot = 7.3e-3;
        fcr3atot = 2.2e-3;
    elseif bodycell == 'penile'
        igg1tot = 2.32;
        igg2tot = 0.0093;
        igg3tot = 0.67;
        igg4tot = 0.00093;
        agtot = 4.2e-7;
        fcr2atot = 2.2e-3;
        fcr3atot = 1.3e-1;
    end
end

%initial conditions for complexes
igg1a_init = 0;
igg2a_init = 0;
igg3a_init = 0;
igg4a_init = 0;

igg1a1_init = 0;
igg1a2_init = 0;
igg1a3_init = 0;
igg1a4_init = 0;

igg2a2_init = 0;
igg2a3_init = 0;
igg2a4_init = 0;

igg3a3_init = 0;
igg3a4_init = 0;

igg4a4_init = 0;

igg1a1r2a_init = 0;
igg1a2r2a_init = 0;
igg1a3r2a_init = 0;
igg1a4r2a_init = 0;

igg2a2r2a_init = 0;
igg2a3r2a_init = 0;
igg2a4r2a_init = 0;

igg3a3r2a_init = 0;
igg3a4r2a_init = 0;

igg4a4r2a_init = 0;

igg1a1r3a_init = 0;
igg1a2r3a_init = 0;
igg1a3r3a_init = 0;
igg1a4r3a_init = 0;

igg2a2r3a_init = 0;
igg2a3r3a_init = 0;
igg2a4r3a_init = 0;

igg3a3r3a_init = 0;
igg3a4r3a_init = 0;

igg4a4r3a_init = 0;

p = [k1af k1ar k2af k2ar k3af k3ar k4af k4ar...
    k1r2f k1r2r k2r2f k2r2r k3r2f k3r2r k4r2f k4r2r...
    k1r3f k1r3r k2r3f k2r3r k3r3f k3r3r k4r3f k4r3r...
    igg1tot igg2tot igg3tot igg4tot agtot fcr2atot fcr3atot];

y0 = [igg1a_init, igg2a_init, igg3a_init, igg4a_init,...
    igg1a1_init, igg1a2_init, igg1a3_init, igg1a4_init,...
    igg2a2_init, igg2a3_init, igg2a4_init,...
    igg3a3_init, igg3a4_init,...
    igg4a4_init...
    igg1a1r2a_init, igg1a2r2a_init, igg1a3r2a_init, igg1a4r2a_init,...
    igg2a2r2a_init, igg2a3r2a_init, igg2a4r2a_init,...
    igg3a3r2a_init, igg3a4r2a_init,...
    igg4a4r2a_init...
    igg1a1r3a_init, igg1a2r3a_init, igg1a3r3a_init, igg1a4r3a_init,...
    igg2a2r3a_init, igg2a3r3a_init, igg2a4r3a_init,...
    igg3a3r3a_init, igg3a4r3a_init,...
    igg4a4r3a_init]';

%Suzie added 2/1
y0 = [igg1a_init, igg2a_init, igg3a_init, igg4a_init,...
    igg1a1_init, igg1a2_init, igg1a3_init, igg1a4_init,...
    igg2a2_init, igg2a3_init, igg2a4_init,...
    igg3a3_init, igg3a4_init,...
    igg4a4_init...
    igg1a1r2a_init, igg1a2r2a_init, igg1a3r2a_init, igg1a4r2a_init,...
    igg2a2r2a_init, igg2a3r2a_init, igg2a4r2a_init,...
    igg3a3r2a_init, igg3a4r2a_init,...
    igg4a4r2a_init...
    igg1a1r3a_init, igg1a2r3a_init, igg1a3r3a_init, igg1a4r3a_init,...
    igg2a2r3a_init, igg2a3r3a_init, igg2a4r3a_init,...
    igg3a3r3a_init, igg3a4r3a_init,...
    igg4a4r3a_init, igg1tot, igg2tot, ...
    igg3tot, igg4tot, agtot, ...
    fcr2atot, fcr3atot]';
% ^^

tspan = [0 100000];

options = odeset('AbsTol',1e-50,'RelTol',1e-10);
% options = odeset('AbsTol',1e-50,'RelTol',1e-10, 'NonNegative',20);
