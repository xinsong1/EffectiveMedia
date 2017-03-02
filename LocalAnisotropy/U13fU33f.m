% f means power of four
function u = U13fU33f(xi) 
K4 = Kint4(xi);
K6 = Kint6(xi);
K8 = Kint8(xi);

u = 3/8*(K4- 2*K6 + K8);
end
