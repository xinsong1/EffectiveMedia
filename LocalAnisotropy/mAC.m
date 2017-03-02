% <AC = <C1111 C3333> 
function AC = mAC(xi, a, b, r, d, e)
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

u13fu33f = U13fU33f(xi);

u13f = U13f(xi);
u33f = U33f(xi);

u13tu33t = U13tU33t(xi);
u13t = U13t(xi);
u33t = U33t(xi);
u13fu33t = U13fU33t(xi);
u13tu33f = U13tU33f(xi);

m  = 2*d+4*e;


AC = a^2 + 4*b^2 + r^2*u13fu33f + m^2*u13tu33t + ...
    4*a*b + a*r*(u13f + u33f) + a*m*(u13t + u33t) + ...
    2*b*r*(u33f + u13f) + 2*b*m*(u13t + u33t) + r*m*(u13fu33t + u13tu33f);
end