% <A^2> = <C1111^2> = <C2222^2>

function AA = mAA(xi, a, b, r, d, e)
% clear all
% close all
% clc
% 

% xi = 0.001;
% 
% Ai = 200;
% Fi = 70;
% Ci = 100;
% Ni = 60;
% Li = 90;
% Bi = Ai-2*Ni;
% 
% a = Ai - 2*Ni;
% b = Ni;
% r = Ai + Ci - 2*Fi - 4*Li;
% d = -Ai + Fi+ 2*Ni;
% e = Li - Ni;

% get needed Tensor U

u13e = U13e(xi);
u13f = U13f(xi);
u13t = U13t(xi);
u13s = U13s(xi);


m = 2*d+4*e;

AA = a^2 + 4*b^2 + r^2*u13e + m^2*u13f + 4*a*b ...
    +2*a*r*u13f + 2*a*m*u13t + ...
    +4*b*r*u13f + 4*b*m*u13t + 2*r*m*u13s;
end