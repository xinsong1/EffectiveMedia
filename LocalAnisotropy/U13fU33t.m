% t means power of two
% f means power of four
function u = U13fU33t(xi) 

K2 = Kint2(xi);
K4 = Kint4(xi);
K6 = Kint6(xi);

u = 3/8*(K2- 2*K4 + K6);
end
