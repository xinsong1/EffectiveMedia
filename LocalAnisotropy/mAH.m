% <AH> = <C1111C2222> 

function AH = mAH(xi, a, b, r, d, e)
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
u13tu23t = U13tU23t(xi);
u13fu23f = U13fU23f(xi);
u13f = U13f(xi);
u23f = U23f(xi);
u13t = U13t(xi);
u23t = U23t(xi);
u13fu23t = U13fU23t(xi);
u13tu23f = u13fu23t;


m = 2*d+4*e;

AH = a^2 + 4*b^2 + r^2*u13fu23f + m^2*u13tu23t + 4*a*b ...
    +a*r*(u13f + u23f) + a*m*(u13t +u23t) + ...
    +2*b*r*(u13f +u23f) + 2*b*m*(u13t+u23t) + r*m*(u13tu23f+u13fu23t);
end