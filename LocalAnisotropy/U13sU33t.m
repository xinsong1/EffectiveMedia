% t means power of two
% s means power of six
function u = U13sU33t(xi) 
K2 = Kint2(xi);
K4 = Kint4(xi);
K6 = Kint6(xi);
K8 = Kint8(xi);

u = 5/16*(K2- 3*K4 + 3*K6 - K8);
end
