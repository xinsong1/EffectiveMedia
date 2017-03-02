function z = tint4(a)
% Closed form of Int[x^4/(a*x^2+1)^2]dx from 0 to 1
% for eta = 1, a =0 very easy
% by Thomas Jordan
% Last modified by Xin Song 03/10/15
if a~= 0;
    z = (-3*atan(sqrt(a))/a^(5/2) + (2*a+3)/(a^3+a^2))/2;
else
    z = 1/5;
end
end
