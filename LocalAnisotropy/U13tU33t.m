% t means power of two
function u = U13tU33t(xi) 

K2 = Kint2(xi);
K4 = Kint4(xi);

u = 1/2*(K2 - K4);
end

