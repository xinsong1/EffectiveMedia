% <AF> = <C1111 C1133> 
function AF = mAF(xi, a, b, r, d, e)
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


u13f = U13f(xi);

u13tu33t = U13tU33t(xi);
u13t = U13t(xi);
u33t = U33t(xi);
u13fu33t = U13fU33t(xi);
u13su33t = U13sU33t(xi);
u13s = U13s(xi);

m  = 2*d+4*e;


AF = a^2 + 2*a*b + a*r*u13f + a*m*u13t + ...
    a*r*u13tu33t + 2*b*r*u13tu33t + r^2 * u13su33t + r*m*u13fu33t + ...
    a*d*(u13t+ u33t) + 2*b*d*(u13t+ u33t) + r*d*(u13s + u13fu33t) + d*m*(u13f + u13tu33t); 
end