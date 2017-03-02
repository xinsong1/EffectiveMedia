% <C>
function C = mC(xi, a, b, r, d, e)




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

u33f = U33f(xi);
u33t = U33t(xi);
C = a + 2*b + r*u33f + (2*d+4*e)*u33t;
end