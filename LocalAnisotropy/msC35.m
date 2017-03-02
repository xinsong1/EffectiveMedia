% <C3313 C1333> = <C1333 C3313>
% msC14 mean squared C14
function sC35 = msC35(xi, r, d, e)
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


sC35 = r^2*u13tu33s + (d^2 + 4*e^2 + 4*d*e)*u13tu33t + (2*r*d + 4*r*e)*u13tu33f; 
end