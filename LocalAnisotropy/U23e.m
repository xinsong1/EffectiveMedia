function u = U23e(xi) 
K0 = Kint0(xi);
K2 = Kint2(xi);
K4 = Kint4(xi);
K6 = Kint6(xi);
K8 = Kint8(xi);
u = 35/128*(K0- 4*K2 + 6*K4 - 4*K6 + K8);
end