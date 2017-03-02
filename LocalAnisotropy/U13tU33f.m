% t means power of two
% f means power of four
function u = U13tU33f(xi) 

K4 = Kint4(xi);
K6 = Kint6(xi);

u = 1/2*(K4 - K6);
end
