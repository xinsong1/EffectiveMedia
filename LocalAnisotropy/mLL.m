% <L^2> = <C1313^2> = <C2323^2> 

function LL = mLL(xi, b, r, e)
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

u13fu33f = U13fU33f(xi);
u13f = U13f(xi);
u33f = U33f(xi);
u13tu33t = U13tU33t(xi);
u13t = U13t(xi);
u33t = U33t(xi);
u13fu33t = U13fU33t(xi);
u13tu33f = U13tU33f(xi);



LL = b^2 + r^2*u13fu33f + e^2*(u13f + u33f + 2*u13tu33t) ...
    +2*b*r*u13tu33t + 2*b*e*(u13t + u33t) + 2*e*r*(u13fu33t + u13tu33f);
end