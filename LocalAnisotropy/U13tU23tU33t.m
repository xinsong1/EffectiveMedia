% t means power of two

function u = U13tU23tU33t(xi) 

K2 = Kint2(xi);
K4 = Kint4(xi);
K6 = Kint6(xi);


u = 1/8*(K2-2*K4+K6);
end
