% <AB> = <C1111 C1122> 
function AB = mAB(xi, a, b, r, d, e)
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

u13tu23t = U13tU23t(xi);
u13t = U13t(xi);
u23t = U23t(xi);
u13fu23t = U13fU23t(xi);

u13su23t = U13sU23t(xi);
u13s = U13s(xi);

m  = 2*d+4*e;


AB = a^2 + 2*a*b + a*r*u13f + a*m*u13t + ...
    a*r*u13tu23t + 2*b*r*u13tu23t + r^2 * u13su23t + r*m*u13fu23t + ...
    a*d*(u13t+ u23t) + 2*b*d*(u13t+ u23t) + r*d*(u13s + u13fu23t) + d*m*(u13f + u13tu23t); 
end