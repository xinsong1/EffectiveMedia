% <C1323 C2313> = <C2313 C1323>
% msC56 mean squared C56
function sC56 = msC56(xi, r, e)
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


sC56 = r^2*u13tu23tu33f + e^2*u13tu23t + 2*r*e*u13tu23tu33t; 
end