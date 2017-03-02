% f means power of four
function u = U23f(xi) 
K0 = Kint0(xi);
K2 = Kint2(xi);
K4 = Kint4(xi);
u = 3/8*(K0-2*K2 + K4);
end