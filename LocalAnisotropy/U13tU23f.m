% t means power of two
% f means power of four
function u = U13tU23f(xi) 
K0 = Kint0(xi);
K2 = Kint2(xi);
K4 = Kint4(xi);
K6 = Kint6(xi);

u = 1/16*(K0- 3*K2 + 3*K4 - K6);
end
