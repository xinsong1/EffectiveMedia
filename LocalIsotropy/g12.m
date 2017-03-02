function z = g12(lambda, mu, eta)
% Thomas H. Jordan 2011/03
% Closed form of spherically averaged STG component g13 = G1133
a = eta^2 - 1;
% if a <= 0, error('eta <= 1'), end
z = (1/(8*mu*(lambda + 2*mu)))*(lambda + mu)*(tint0(a) -2*tint2(a)+ tint4(a));