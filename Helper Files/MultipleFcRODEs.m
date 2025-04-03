function [dydt, freeSpec] = MultipleFcRODEs(t,y,p)
k1af = p(1);
k1ar = p(2);
k2af = p(3);
k2ar = p(4);
k3af = p(5);
k3ar = p(6);
k4af = p(7);
k4ar = p(8);

k1r2f = p(9);
k1r2r = p(10);
k2r2f = p(11);
k2r2r = p(12);
k3r2f = p(13);
k3r2r = p(14);
k4r2f = p(15);
k4r2r = p(16);

k1r3f = p(17);
k1r3r = p(18);
k2r3f = p(19);
k2r3r = p(20);
k3r3f = p(21);
k3r3r = p(22);
k4r3f = p(23);
k4r3r = p(24);

igg1tot = p(25);
igg2tot = p(26);
igg3tot = p(27);
igg4tot = p(28);
agtot = p(29);
fcr2atot = p(30);
fcr3atot = p(31);

igg1a = y(1);
igg2a = y(2);
igg3a = y(3);
igg4a = y(4);

igg1a1 = y(5);
igg1a2 = y(6);
igg1a3 = y(7);
igg1a4 = y(8);
igg2a2 = y(9);
igg2a3 = y(10);
igg2a4 = y(11);
igg3a3 = y(12);
igg3a4 = y(13);
igg4a4 = y(14);

igg1a1r2a = y(15);
igg1a2r2a = y(16);
igg1a3r2a = y(17);
igg1a4r2a = y(18);
igg2a2r2a = y(19);
igg2a3r2a = y(20);
igg2a4r2a = y(21);
igg3a3r2a = y(22);
igg3a4r2a = y(23);
igg4a4r2a = y(24);

igg1a1r3a = y(25);
igg1a2r3a = y(26);
igg1a3r3a = y(27);
igg1a4r3a = y(28);
igg2a2r3a = y(29);
igg2a3r3a = y(30);
igg2a4r3a = y(31);
igg3a3r3a = y(32);
igg3a4r3a = y(33);
igg4a4r3a = y(34);

ag = agtot - igg1a - igg2a - igg3a - igg4a...
    - igg1a1 - igg1a2-  igg1a3 - igg1a4...
    - igg2a2 - igg2a3 - igg2a4...
    - igg3a3 - igg3a4...
    - igg4a4...
    - igg1a1r2a - igg1a2r2a - igg1a3r2a - igg1a4r2a...
    - igg2a2r2a - igg2a3r2a - igg2a4r2a...
    - igg3a3r2a - igg3a4r2a...
    - igg4a4r2a...
    - igg1a1r3a - igg1a2r3a - igg1a3r3a - igg1a4r3a...
    - igg2a2r3a - igg2a3r3a - igg2a4r3a...
    - igg3a3r3a - igg3a4r3a...
    - igg4a4r3a;

igg1 = igg1tot - igg1a - 2*igg1a1 - igg1a2 - igg1a3 - igg1a4...
    - 2*igg1a1r2a - igg1a2r2a - igg1a3r2a - igg1a4r2a...
    - 2*igg1a1r3a - igg1a2r3a - igg1a3r3a - igg1a4r3a;

igg2 = igg2tot - igg2a - igg1a2 - 2*igg2a2 - igg2a3 - igg2a4...
    - igg1a2r2a - 2*igg2a2r2a - igg2a3r2a - igg2a4r2a...
    - igg1a2r3a - 2*igg2a2r3a - igg2a3r3a - igg2a4r3a;

igg3 = igg3tot - igg3a - igg1a3 - igg2a3 - 2*igg3a3 - igg3a4...
    - igg1a3r2a - igg2a3r2a - 2*igg3a3r2a - igg3a4r2a...
    - igg1a3r3a - igg2a3r3a - 2*igg3a3r3a - igg3a4r3a;

igg4 = igg4tot - igg4a - igg1a4 - igg2a4 - igg3a4 - 2*igg4a4...
    - igg1a4r2a - igg2a4r2a - igg3a4r2a - 2*igg4a4r2a...
    - igg1a4r3a - igg2a4r3a - igg3a4r3a - 2*igg4a4r3a;

fcr2a = fcr2atot...
    - igg1a1r2a - igg1a2r2a - igg1a3r2a - igg1a4r2a...
    - igg2a2r2a - igg2a3r2a - igg2a4r2a...
    - igg3a3r2a - igg3a4r2a...
    - igg4a4r2a;

fcr3a =fcr3atot...
    - igg1a1r3a - igg1a2r3a - igg1a3r3a - igg1a4r3a...
    - igg2a2r3a - igg2a3r3a - igg2a4r3a...
    - igg3a3r3a - igg3a4r3a...
    - igg4a4r3a;

%ODEs
d1a = 2*k1af*igg1*ag - k1ar*igg1a...
    - k1af*igg1a*igg1 + 2*k1ar*igg1a1...
    - k2af*igg1a*igg2 + k2ar*igg1a2...
    - k3af*igg1a*igg3 + k3ar*igg1a3...
    - k4af*igg1a*igg4 + k4ar*igg1a4;

d2a = 2*k2af*igg2*ag - k2ar*igg2a...
    - k1af*igg2a*igg1 + k1ar*igg1a2...
    - k2af*igg2a*igg2 + 2*k2ar*igg2a2...
    - k3af*igg2a*igg3 + k3ar*igg2a3...
    - k4af*igg2a*igg4 + k4ar*igg2a4;

d3a = 2*k3af*igg3*ag - k3ar*igg3a...
    - k1af*igg3a*igg1 + k1ar*igg1a3...
    - k2af*igg3a*igg2 + k2ar*igg2a3...
    - k3af*igg3a*igg3 + 2*k3ar*igg3a3...
    - k4af*igg3a*igg4 + k4ar*igg3a4;

d4a = 2*k4af*igg4*ag - k4ar*igg4a...
    - k1af*igg4a*igg1 + k1ar*igg1a4...
    - k2af*igg4a*igg2 + k2ar*igg2a4...
    - k3af*igg4a*igg3 + k3ar*igg3a4...
    - k4af*igg4a*igg4 + 2*k4ar*igg4a4;

d1a1 = k1af*igg1a*igg1 - 2*k1ar*igg1a1...
    - k1r2f*igg1a1*fcr2a + k1r2r*igg1a1r2a...
    - k1r3f*igg1a1*fcr3a + k1r3r*igg1a1r3a;

d1a2 = k2af*igg1a*igg2 - k2ar*igg1a2...
    + k1af*igg2a*igg1 - k1ar*igg1a2...
    - mean([k1r2f k2r2f])*igg1a2*fcr2a + mean([k1r2r k2r2r])*igg1a2r2a...
    - mean([k1r3f k2r3f])*igg1a2*fcr3a + mean([k1r3r k2r3r])*igg1a2r3a;

d1a3 = k3af*igg1a*igg3 - k3ar*igg1a3...
    + k1af*igg3a*igg1 - k1ar*igg1a3...
    - mean([k1r2f k3r2f])*igg1a3*fcr2a + mean([k1r2r k3r2r])*igg1a3r2a...
    - mean([k1r3f k3r3f])*igg1a3*fcr3a + mean([k1r3r k3r3r])*igg1a3r3a;

d1a4 = k4af*igg1a*igg4 - k4ar*igg1a4...
    + k1af*igg4a*igg1 - k1ar*igg1a4...
    - mean([k1r2f k4r2f])*igg1a4*fcr2a + mean([k1r2r k4r2r])*igg1a4r2a...
    - mean([k1r3f k4r3f])*igg1a4*fcr3a + mean([k1r3r k4r3r])*igg1a4r3a;

d2a2 = k2af*igg2a*igg2 - 2*k2ar*igg2a2...
    - k2r2f*igg2a2*fcr2a + k2r2r*igg2a2r2a...
    - k2r3f*igg2a2*fcr3a + k2r3r*igg2a2r3a;

d2a3 = k3af*igg2a*igg3 - k3ar*igg2a3...
    + k2af*igg3a*igg2 - k2ar*igg2a3...
    - mean([k2r2f k3r2f])*igg2a3*fcr2a + mean([k2r2r k3r2r])*igg2a3r2a...
    - mean([k2r3f k3r3f])*igg2a3*fcr3a + mean([k2r3r k3r3r])*igg2a3r3a;

d2a4 = k4af*igg2a*igg4 - k4ar*igg2a4...
    + k2af*igg4a*igg2 - k2ar*igg2a4...
    - mean([k2r2f k4r2f])*igg2a4*fcr2a + mean([k2r2r k4r2r])*igg2a4r2a...
    - mean([k2r3f k4r3f])*igg2a4*fcr3a + mean([k2r3r k4r3r])*igg2a4r3a;

d3a3 = k3af*igg3a*igg3 - 2*k3ar*igg3a3...
    - k3r2f*igg3a3*fcr2a + k3r2r*igg3a3r2a...
    - k3r3f*igg3a3*fcr3a + k3r3r*igg3a3r3a;

d3a4 = k4af*igg3a*igg4 - k4ar*igg3a4...
    + k3af*igg4a*igg3 - k3ar*igg3a4...
    - mean([k3r2f k4r2f])*igg3a4*fcr2a + mean([k3r2r k4r2r])*igg3a4r2a...
    - mean([k3r3f k4r3f])*igg3a4*fcr3a + mean([k3r3r k4r3r])*igg3a4r3a;

d4a4 = k4af*igg4a*igg4 - 2*k4ar*igg4a4...
    - k4r2f*igg4a4*fcr2a + k4r2r*igg4a4r2a...
    - k4r3f*igg4a4*fcr3a + k4r3r*igg4a4r3a;

d1a1r2a = k1r2f*igg1a1*fcr2a - k1r2r*igg1a1r2a;

d1a2r2a = mean([k1r2f k2r2f])*igg1a2*fcr2a - mean([k1r2r k2r2r])*igg1a2r2a;

d1a3r2a = mean([k1r2f k3r2f])*igg1a3*fcr2a - mean([k1r2r k3r2r])*igg1a3r2a;

d1a4r2a = mean([k1r2f k4r2f])*igg1a4*fcr2a - mean([k1r2r k4r2r])*igg1a4r2a;

d2a2r2a = k2r2f*igg2a2*fcr2a - k2r2r*igg2a2r2a;

d2a3r2a = mean([k2r2f k3r2f])*igg2a3*fcr2a - mean([k2r2r k3r2r])*igg2a3r2a;

d2a4r2a = mean([k2r2f k4r2f])*igg2a4*fcr2a - mean([k2r2r k4r2r])*igg2a4r2a;

d3a3r2a = k3r2f*igg3a3*fcr2a - k3r2r*igg3a3r2a;

d3a4r2a = mean([k3r2f k4r2f])*igg3a4*fcr2a - mean([k3r2r k4r2r])*igg3a4r2a;

d4a4r2a = k4r2f*igg4a4*fcr2a - k4r2r*igg4a4r2a;

d1a1r3a = k1r3f*igg1a1*fcr3a - k1r3r*igg1a1r3a;

d1a2r3a = mean([k1r3f k2r3f])*igg1a2*fcr3a - mean([k1r3r k2r3r])*igg1a2r3a;

d1a3r3a = mean([k1r3f k3r3f])*igg1a3*fcr3a - mean([k1r3r k3r3r])*igg1a3r3a;

d1a4r3a = mean([k1r3f k4r3f])*igg1a4*fcr3a - mean([k1r3r k4r3r])*igg1a4r3a;

d2a2r3a = k2r3f*igg2a2*fcr3a - k2r3r*igg2a2r3a;

d2a3r3a = mean([k2r3f k3r3f])*igg2a3*fcr3a - mean([k2r3r k3r3r])*igg2a3r3a;

d2a4r3a = mean([k2r3f k4r3f])*igg2a4*fcr3a - mean([k2r3r k4r3r])*igg2a4r3a;

d3a3r3a = k3r3f*igg3a3*fcr3a - k3r3r*igg3a3r3a;

d3a4r3a = mean([k3r3f k4r3f])*igg3a4*fcr3a - mean([k3r3r k4r3r])*igg3a4r3a;

d4a4r3a = k4r3f*igg4a4*fcr3a - k4r3r*igg4a4r3a;

dydt = [d1a d2a d3a d4a...
    d1a1 d1a2 d1a3 d1a4...
    d2a2 d2a3 d2a4...
    d3a3 d3a4...
    d4a4...
    d1a1r2a d1a2r2a d1a3r2a d1a4r2a...
    d2a2r2a d2a3r2a d2a4r2a...
    d3a3r2a d3a4r2a...
    d4a4r2a...
    d1a1r3a d1a2r3a d1a3r3a d1a4r3a...
    d2a2r3a d2a3r3a d2a4r3a...
    d3a3r3a d3a4r3a...
    d4a4r3a]';

freeSpec = [igg1 igg2 igg3 igg4 ag fcr2a fcr3a];

% added free species
dydt = [dydt; freeSpec'];

% disp('Still running')