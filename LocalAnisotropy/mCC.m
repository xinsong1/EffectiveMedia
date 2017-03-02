% <C^2> = <C3333^2> 

function CC = mCC(xi, a, b, r, d, e)
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

u33e = U33e(xi);
u33f = U33f(xi);
u33t = U33t(xi);
u33s = U33s(xi);


m = 2*d+4*e;

CC = a^2 + 4*b^2 + r^2*u33e + m^2*u33f + 4*a*b ...
    +2*a*r*u33f + 2*a*m*u33t + ...
    +4*b*r*u33f + 4*b*m*u33t + 2*r*m*u33s;
end