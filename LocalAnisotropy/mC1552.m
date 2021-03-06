% <C1113 C1322> = <C2213 C1311>
function C1552 = mC1552(xi, r, d, e)
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
u13fu33t = U13fU33t(xi);
u13tu33t = U13tU33t(xi);


C1552 = r^2*u13fu23tu33t + (d*r+ 2*e*r)*u13tu23tu33t + r*d*u13fu33t + (d^2+2*e*d)* u13tu33t; 
end