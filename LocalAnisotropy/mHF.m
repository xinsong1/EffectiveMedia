% <HF> = <C2222 C1133> 
function HF = mHF(xi, a, b, r, d, e)
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


u23f = U23f(xi);
u13tu23t = U13tU23t(xi);
u13tu33t = U13tU33t(xi);
u13t = U13t(xi);
u23t = U23t(xi);
u33t = U33t(xi);
u13fu33t = U13fU33t(xi);
u23fu33t = u13fu33t;
u13fu23t = U13fU23t(xi);
u13tu23f = u13fu23t;
u13fu23tu33t = U13fU23tU33t(xi);
u13tu23fu33t = u13fu23tu33t;
u13tu23tu33t = U13tU23tU33t(xi);


m  = 2*d+4*e;


HF = a^2 + 2*a*b + a*r*u23f + a*m*u23t + ...
    a*r*u13tu33t + 2*b*r*u13tu33t + r^2 * u13tu23fu33t + r*m*u13tu23tu33t + ...
    a*d*(u13t+ u33t) + 2*b*d*(u13t+ u33t) + r*d*(u13tu23f + u23fu33t) + d*m*(u13tu23t + u13tu33t); 
end