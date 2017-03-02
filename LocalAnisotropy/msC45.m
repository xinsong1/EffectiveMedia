% <C1312 C1213> = <C1213 C1312>
% msC45 mean squared C45
function sC45 = msC45(xi, r, e)
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
u13tu33t = U13tU33t(xi);
u23tu33t = u13tu33t;
u13tu23tu33t = U13tU23tU33t(xi);


sC45 = r^2*u13fu23tu33t + e^2*u23tu33t + 2*r*e*u13tu23tu33t; 
end