% t means power of two
function u = U13tU23t(xi) 
K0 = Kint0(xi);
K2 = Kint2(xi);
K4 = Kint4(xi);

u = 1/8*(K0- 2*K2 + K4);
end

