% <C1113 C1311> = <C1311 C1113>
% msC15 mean squared C15
function sC15 = msC15(xi, r, d, e)
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
u13su33t = U13sU33t(xi);
u13tu33t = U13tU33t(xi);
u13fu33t = U13fU33t(xi);


sC15 = r^2*u13su33t + (d^2 + 4*e^2 + 4*d*e)*u13tu33t + (2*r*d + 4*r*e)*u13fu33t; 
end