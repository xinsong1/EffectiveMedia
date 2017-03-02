% <C1223 C2312> = <C2312 C1223>
% msC46 mean squared C46
function sC46 = msC46(xi, r, e)
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
u13tu23fu33t = u13fu23tu33t;
u13tu33t = U13tU33t(xi);
u13tu23tu33t = U13tU23tU33t(xi);


sC46 = r^2*u13tu23fu33t + e^2*u13tu33t + 2*r*e*u13tu23tu33t; 
end