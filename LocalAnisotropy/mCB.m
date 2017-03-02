% <CB> = <C3333 C1122> 
function CB = mCB(xi, a, b, r, d, e)
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

u33t = U33t(xi);
u33f = U33f(xi);

u13tu23t = U13tU23t(xi);
u13tu33t = U13tU33t(xi);
u23tu33t = u13tu33t;
u13t = U13t(xi);

u23t = U23t(xi);
u13tu23tu33t = U13tU23tU33t(xi);
u13tu23tu33f = U13tU23tU33f(xi);
u13tu33f = U13tU33f(xi);
u23tu33f = u13tu33f;
m  = 2*d+4*e;


CB = a^2 + 2*a*b + a*r*u33f + a*m*u33t + ...
    a*r*u13tu23t + 2*b*r*u13tu23t + r^2 * u13tu23tu33f + r*m*u13tu23tu33t + ...
    a*d*(u13t+ u23t) + 2*b*d*(u13t+ u23t) + r*d*(u13tu33f + u23tu33f) + d*m*(u13tu33t + u23tu33t); 
end