function z = tint2(a)
% Closed form of Int[x^2/(a*x^2+1)^2]dx from 0 to 1
% for eta = 1, a = 0 very easy
% by Thomas Jordan
% Last modified by Xin Song 03/10/15
if a ~= 0 
    z = (atan(sqrt(a))/a^(3/2) - 1/(a^2+a))/2;
else
    z = 1/3;
end
end
