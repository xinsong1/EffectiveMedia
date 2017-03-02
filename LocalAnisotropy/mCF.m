% <CF> = <C3333 C1133> 
function CF = mCF(xi, a, b, r, d, e)
% clear all
% close all
% clc


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


u33f = U33f(xi);

u13tu33t = U13tU33t(xi);
u13t = U13t(xi);
u33t = U33t(xi);
u13tu33f = U13tU33f(xi);
u13tu33s = U13tU33s(xi);
u33s = U33s(xi);

m  = 2*d+4*e;


CF = a^2 + 2*a*b + a*r*u33f + a*m*u33t + ...
    a*r*u13tu33t + 2*b*r*u13tu33t + r^2 * u13tu33s + r*m*u13tu33f + ...
    a*d*(u13t+ u33t) + 2*b*d*(u13t+ u33t) + r*d*(u33s + u13tu33f) + d*m*(u33f + u13tu33t); 
end