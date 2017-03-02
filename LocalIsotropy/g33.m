function z = g33(lambda, mu, eta)
% Thomas H. Jordan 2011/03
% Closed form of spherically averaged STG component g33 = G3333
a = eta^2 - 1;
% if a <= 0, error('eta <= 1'), end
z = -(1/(mu*(lambda + 2*mu)))*(a + 1)*((a + 1)*mu*tint4(a) + (lambda + 2*mu)*(tint2(a) - tint4(a)));