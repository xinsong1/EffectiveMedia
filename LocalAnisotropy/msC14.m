% <C1112 C1211> = <C1211 C1112>
% msC14 mean squared C14
function sC14 = msC14(xi, r, d, e)
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
u13su23t = U13sU23t(xi);
u13tu23t = U13tU23t(xi);
u13fu23t = U13fU23t(xi);


sC14 = r^2*u13su23t + (d^2 + 4*e^2 + 4*d*e)*u13tu23t + (2*r*d + 4*r*e)*u13fu23t; 
end