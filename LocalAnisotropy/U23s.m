% s means power of six
function u = U23s(xi) 
K0 = Kint0(xi);
K2 = Kint2(xi);
K4 = Kint4(xi);
K6 = Kint6(xi);
u = 5/16*(K0-3*K2 + 3*K4 - K6);
end