function z = Kint2(xi)
% Calculate Eq 53 from Jordan 2015
% n = 2
% Xin Song 03/30/2015
if xi > 1
    a = xi^2-1;
    z = 1/a*(xi/sqrt(a)*asinh(sqrt(a))-1);
elseif xi ==1
    z = 1/3;
elseif xi > 0 && xi < 1
    a = 1-xi^2;
    z = 1/a*(1-xi*acos(xi)/a^(1/2));
elseif xi == 0
    z = 1;
end

end