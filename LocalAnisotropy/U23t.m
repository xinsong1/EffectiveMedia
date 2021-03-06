% Eq C8b
% t means power of two
function u = U23t(xi) 
K0 = Kint0(xi);
K2 = Kint2(xi);
u = 1/2*(K0-K2);
end