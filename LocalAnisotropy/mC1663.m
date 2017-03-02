% <C1123 C2333> = <C3323 C2311>
function C1663 = mC1663(xi, r, d, e)
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
u13tu23tu33f = U13tU23tU33f(xi);
u13tu23tu33t = U13tU23tU33t(xi);

u13tu33f = U13tU33f(xi);
u23tu33f = u13tu33f;

u13tu33t = U13tU33t(xi);
u23tu33t = u13tu33t;


C1663 = r^2*u13tu23tu33f + (d*r+ 2*e*r)*u13tu23tu33t + r*d*u23tu33f + (d^2+2*e*d)* u23tu33t; 
end