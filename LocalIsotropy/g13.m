function z = g13(lambda, mu, eta)
% Thomas H. Jordan 2011/03
% Closed form of spherically averaged STG component g13 = G1133
a = eta^2 - 1;
% if a <= 0, error('eta <= 1'), end
z = (1/(2*mu*(lambda + 2*mu)))*(lambda + mu)*(a + 1)*(tint2(a) - tint4(a));