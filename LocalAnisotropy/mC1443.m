% <C1112 C1233> = <C3312 C1211>
function C1443 = mC1443(xi, r, d, e)
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
u13tu23tu33t = U13tU23tU33t(xi);
u13fu23t = U13fU23t(xi);
u13tu23t = U13tU23t(xi);


C1443 = r^2*u13fu23tu33t + (d*r+ 2*e*r)*u13tu23tu33t + r*d*u13fu23t + (d^2+2*e*d)* u13tu23t; 
end