% t means power of two
% s means power of six
function u = U13sU23t(xi) 
K0 = Kint0(xi);
K2 = Kint2(xi);
K4 = Kint4(xi);
K6 = Kint6(xi);
K8 = Kint8(xi);

u = 5/128*(K0- 4*K2 + 6*K4 - 4*K6 + K8);
end
