% <A>
function A = mA(xi, a, b, r, d, e)




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

u13f = U13f(xi);
u13t = U13t(xi);
A = a + 2*b + r*u13f + (2*d+4*e)*u13t;
end