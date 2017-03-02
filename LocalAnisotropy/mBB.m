% <B^2> = <C1122^2> = <C2211^2>

function BB = mBB(xi, a, r, d)
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

u13fu23f = U13fU23f(xi);

u13f = U13f(xi);
u23f = U23f(xi);
u13tu23t = U13tU23t(xi);
u13t = U13t(xi);
u23t = U23t(xi);
u13fu23t = U13fU23t(xi);
u13tu23f = U13tU23f(xi);



BB = a^2 + r^2*u13fu23f + d^2*(u13f + u23f + 2*u13tu23t) ...
    +2*a*r*u13tu23t + 2*a*d*(u13t + u23t) + 2*r*d*(u13fu23t + u13tu23f);
  
end