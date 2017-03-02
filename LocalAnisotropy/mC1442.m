% <C1112 C1222> = <C2212 C1211>
function C1442 = mC1442(xi, r, d, e)
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
u13fu23f = U13fU23f(xi);
u13tu23t = U13tU23t(xi);
u13fu23t = U13fU23t(xi);
u13tu23f = U13tU23f(xi);


C1442 = r^2*u13fu23f + (d^2 + 4*e^2 + 4*d*e)*u13tu23t + (r*d + 2*r*e)*(u13fu23t + u13tu23f); 
end