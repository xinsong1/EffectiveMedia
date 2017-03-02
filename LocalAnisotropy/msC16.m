% <C1123 C2311> = <C2311 C1123>
% msC16 mean squared C44
function sC16 = msC16(xi, r, d)
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
u13fu23tu33t = U13fU23tU33t(xi);
u13tu23t = U13tU23t(xi);
u13tu23tu33t = U13tU23tU33t(xi);


sC16 = r^2*u13fu23tu33t + d^2 *u13tu23t + 2*r*d*u13tu23tu33t; 
end