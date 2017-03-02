function z = tint0(a)
% Closed form of Int[1/(a*x^2+1)^2]dx from 0 to 1
% for eta = 1, a =0 very easy
% by Thomas Jordan
% Last modified by Xin Song 03/10/15
if a ~= 0
    z = (atan(sqrt(a))/sqrt(a) + 1/(a+1))/2;
else
    z = 1;
end
end
