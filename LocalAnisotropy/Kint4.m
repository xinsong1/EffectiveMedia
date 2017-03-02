function z = Kint4(xi, factor)
% Calculate Eq 53 from Jordan 2015
% n = 4
% for xi > 1, two cases:
%               one is the directly integral of Eq 54
%               second one is from Eq C10
% same value
% Xin Song 03/30/2015
% Modified by Xin Song 03/31/2015

if nargin < 2
    factor = 1;
end
if xi > 1
    a = xi.^2-1;
    switch factor

        case 1
            z = 1/2/a*(1-3*Kint2(xi));

        case 2
            z = 1i*(sqrt(a)*(2+xi^2) - 3*xi*asinh(sqrt(a)))/2/sqrt(-a)^5;
    end
elseif xi == 1
    z = 1/5;
elseif xi < 1 && xi > 0
    a = 1-xi^2;
    z = xi*(2/xi+xi+3i*acosh(xi)/sqrt(a))/2/(-a)^2;
elseif xi == 0
    z = 1;
end
end