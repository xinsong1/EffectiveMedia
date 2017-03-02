% <C1113 C1333> = <C3313 C1311>
function C1553 = mC1553(xi, r, d, e)
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
u13fu33f = U13fU33f(xi);
u13tu33t = U13tU33t(xi);
u13fu33t = U13fU33t(xi);
u13tu33f = U13tU33f(xi);


C1553 = r^2*u13fu33f + (d^2 + 4*e^2 + 4*d*e)*u13tu33t + (r*d + 2*r*e)*(u13fu33t + u13tu33f); 
end