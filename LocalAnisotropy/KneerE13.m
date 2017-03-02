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
function z = KneerE13(eta, A, C, F, L)
fun = @(x) (F+L)*eta^2*x.^2.*(1-x.^2)./(2*(C*L*eta^4.*x.^4 + (A*C-F^2-2*F*L)*eta^2.*x.^2.*(1-x.^2)+ A*L.*(1-x.^2).^2));
z = integral(fun, 0, 1);
end

% fun = @(t) (80* (cos(t)).^2 .* (sin(t)).^2)./ ...
%     (9000 * (cos(t)).^4 + 2500* (cos(t)).^2 .* (sin(t)).^2 + 18000* (sin(t)).^4) .* sin(t);
% b = 0.5*integral(fun, 0, pi);