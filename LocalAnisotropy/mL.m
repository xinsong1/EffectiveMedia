% <L>
function L = mL(xi, b, r, e)




% Ai = 200;
% Fi = 50;
% Ci = 100;
% Ni = 60;
% Li = 90;
% Bi = Ai-2*Ni;
% 
% xi = 0.01;
% 
% 
% a = Ai - 2*Ni;
% b = Ni;
% r = Ai + Ci - 2*Fi - 4*Li;
% d = -Ai + Fi+ 2*Ni;
% e = Li - Ni;


u13t = U13t(xi);
u33t = U33t(xi);
u13tu33t = U13tU33t(xi);
L = b + r*u13tu33t + e*(u13t + u33t);
end