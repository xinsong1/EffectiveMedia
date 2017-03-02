% <BF> = <C1122 C1133> 
function BF = mBF(xi, a, r, d)
% clear all
% close all
% clc
% 
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


u13f = U13f(xi);

u13tu33t = U13tU33t(xi);
u23tu33t = u13tu33t;

u13tu23t = U13tU23t(xi);
u13t = U13t(xi);
u23t = u13t;
u33t = U33t(xi);

u13fu23tu33t = U13fU23tU33t(xi);
u13tu23tu33t = U13tU23tU33t(xi);

u13fu33t = U13fU33t(xi);
u13fu23t = U13fU23t(xi);


BF = a^2 + r^2*u13fu23tu33t + d^2*(u13f + u13tu23t + u13tu33t + u23tu33t) + ...
    a*r*(u13tu33t + u13tu23t) + a*d*(u23t + u33t + 2*u13t) + ...
    r*d*(u13fu33t+ 2*u13tu23tu33t + u13fu23t);
end