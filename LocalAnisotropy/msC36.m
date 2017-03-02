% <C1112 C1211> = <C1211 C1112>
% msC36 mean squared C36
function sC36 = msC36(xi, r, d, e)
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
u13tu33s = U13tU33s(xi);
u13tu33t = U13tU33t(xi);
u13tu33f = U13tU33f(xi);

u23tu33s = u13tu33s;
u23tu33t = u13tu33t;
u23tu33f = u13tu33f;

sC36 = r^2*u23tu33s + (d^2 + 4*e^2 + 4*d*e)*u23tu33t + (2*r*d + 4*r*e)*u23tu33f; 
end