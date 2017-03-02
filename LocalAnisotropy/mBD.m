% <BF> = <C1122 C1133> 
function BD = mBD(xi, a, r, d)
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


u23f = U23f(xi);

u13tu33t = U13tU33t(xi);
u23tu33t = u13tu33t;

u13tu23t = U13tU23t(xi);
u13t = U13t(xi);
u23t = u13t;
u33t = U33t(xi);

u13fu23tu33t = U13fU23tU33t(xi);
u13tu23fu33t = u13fu23tu33t;
u13tu23tu33t = U13tU23tU33t(xi);

u13fu33t = U13fU33t(xi);
u23fu33t = u13fu33t;
u13fu23t = U13fU23t(xi);
u13tu23f = u13fu23t;


BD = a^2 + r^2*u13tu23fu33t + d^2*(u23f + u13tu23t + u13tu33t + u23tu33t) + ...
    a*r*(u23tu33t + u13tu23t) + a*d*(2*u23t + u33t + u13t) + ...
    r*d*(u23fu33t+ 2*u13tu23tu33t + u13tu23f);
end