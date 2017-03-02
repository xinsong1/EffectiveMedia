% <C3312 C1233> = <C1233 C3312>
% msC34 mean squared C34
function sC34 = msC34(xi, r, d)
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
u13tu23t = U13tU23t(xi);
u13tu23tu33t = U13tU23tU33t(xi);


sC34 = r^2*u13tu23tu33f + d^2 *u13tu23t + 2*r*d*u13tu23tu33t; 
end