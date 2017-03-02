function z = Kint6(xi, factor)
% Calculate Eq 53 from Jordan 2015
% n = 6
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
            z = 1/4/a*(1-5*Kint4(xi, 1));
         case 2
             z = -1i*(sqrt(a)*(-8-9*xi^2+2*xi^4)+15*xi*asinh(sqrt(a)))/8/sqrt(-a)^7;
     end
 
 elseif xi == 1
     z = 1/7;
 
 elseif xi > 0 && xi < 1
     a = 1-xi^2;
     z = -xi*(8/xi + 9*xi - 2*xi^3 + 15i*acosh(xi)/sqrt(a))/8/(-a)^3;
     
 elseif xi == 0
     z = 1;
 end
 
 
    
end