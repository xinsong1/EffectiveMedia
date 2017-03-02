% clear; close all; clc;
% 
% A = 200;
% F = 70;
% C = 100;
% N = 60;
% L = 90;
% B = A-2*N;
% 
% eta = 1;
function z = KneerE11(eta, A, C, F, L, N)
fun = @(x)-(1-x.^2).*(3*(C*eta^2*x.^2 + L*(1-x.^2))+ ...
    (C*L*eta^4*x.^4 + (A*C - F^2 - 2*F*L)*eta^2.*x.^2.*(1-x.^2) + A*L*(1-x.^2).^2)./(L*eta^2*x.^2 + N*(1-x.^2)))./ ...
    (8*(C*L*eta^4.*x.^4 + (A*C-F^2-2*F*L)*eta^2.*x.^2.*(1-x.^2)+ A*L.*(1-x.^2).^2));
z = integral(fun, 0, 1);
end
