function z = Kint8(xi, factor)
% Calculate Eq 53 from Jordan 2015
% n = 8
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
        a = xi^2-1;
        switch factor
            
            case 1
                z = 1/6/a*(1-7*Kint6(xi, 1));

            case 2
                z = 1i*(sqrt(a)*(48+87*xi^2-38*xi^4+8*xi^6)-105*xi*asinh(sqrt(a)))/48/sqrt(-a)^9;
        end

    elseif xi == 1
        z = 1/9;

    elseif xi < 1 && xi > 0
        a = 1-xi^2;
        z = xi*(48/xi + 87*xi - 38*xi^3 + 8*xi^5 + 105i*acosh(xi)/sqrt(a))/48/(-a)^4;
    elseif xi == 0
        z = 1;
    end
end