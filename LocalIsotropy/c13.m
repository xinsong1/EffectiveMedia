function F = c13(lambda, mu, eta, epsl, epsm, alpha, factor, qratio)
% Closed form of effective medium component c13 = C1133 = F
% Elastic parameter calculation:
            % Factor = 1----------Jordan's EM 
            % Factor = 2----------Backus Approximate Theory
            % Factor = 3----------Backus Exact Theory (Gamma) lambda and mu
            % are i.i.d
            % Factor = 4----------Backus Exact Theory (Uniform) lambda and
            % mu are i.i.d
            % Factor = 5----------Backus exact Theory (Gamma, mu =
            % fact*lambda
            % Factor = 6----------Backus exact Theory (Uniform, mu =
            % fact*lambda
            % Factor = 7----------Backus exact Theory (Uniform, qratio involved)
% Assumes dlambda/lambda = dmu/mu = epsilon

% By Thomas Jordan 2011/03 for case 1
% modified by Xin Song 10/27/2014 11/06/2014 2/21/2016 2/22/2016 2/23/2016


if lambda <= 0, error('lambda <= 1'), end
if mu <= 0, error('mu <= 1'), end
if epsl < 0, error('epsl <= 0'), end
if epsm < 0, error('epsm <= 0'), end
% if eta <= 1, error('eta <= 1'), end
dlambda = epsl*lambda;
dmu = epsm*mu;
dlm = alpha*dmu*dlambda;

switch factor
    case 1 % Jordan Effective Medium Theory

        
        F = lambda...
            + 4*(dlambda^2 + dlm)*g11(lambda,mu,eta)...
            + 2*(2*dlambda^2 + 3*dlm + 2*dmu^2)*g13(lambda,mu,eta)...
            + (dlambda^2 + 2*dlm)*g33(lambda,mu,eta)...
            - 4*(dlambda^2 + dlm)*g44(lambda,mu,eta);
    case 2 % Backus 2nd Order Approxiamte Theory
         F = lambda*...
            (1 - (dlambda^2+2*dlm)/(lambda*(lambda+2*mu)));
        
    case 3 % Backus Exact Theory (Gamma) lambda and mu are i.i.d
        % for mu
        ax = mu^2/dmu^2;
        invbx = mu/dmu^2;
        % for lambda
        ay = lambda^2/dlambda^2;
        invby = lambda/dlambda^2;
    
        % integral over gamma distribution
        fun1 = @(x,y) 1./(y+2*x).*1./gamma(ax).*invbx.^ax.*x.^(ax-1).*exp(-invbx.*x) ...
            .*1./gamma(ay).*invby.^ay.*y.^(ay-1).*exp(-invby.*y);
        F1 = integral2(fun1, 0, Inf, 0, Inf);
        F1 = 1/F1;
        
        fun2 = @(x,y) y./(y+2*x).*1./gamma(ax).*invbx.^ax.*x.^(ax-1).*exp(-invbx.*x) ...
            .*1./gamma(ay).*invby.^ay.*y.^(ay-1).*exp(-invby.*y);
        
        F2 = integral2(fun2, 0, Inf, 0, Inf);
        F = F1*F2;
        
        
    case 4 % Backus Exact Theory (Uniform) lambda and mu are i.i.d

        lowerx = (2*mu-sqrt(12)*dmu)/2; upperx = (2*mu+sqrt(12)*dmu)/2; %for mu
        lowery = (2*lambda-sqrt(12)*dlambda)/2; uppery = (2*lambda+sqrt(12)*dlambda)/2; % for lambda
        
        % integral over gamma distribution
        fun1 = @(x,y) 1./(y+2*x).*(upperx-lowerx)^-1.*(uppery-lowery)^-1;
        F1 = integral2(fun1, lowerx, upperx, lowery, uppery);
        F1 = 1/F1;
        
        fun2 = @(x,y) y./(y+2*x).*(upperx-lowerx)^-1.*(uppery-lowery)^-1;
        
        F2 = integral2(fun2, lowerx, upperx, lowery, uppery);
        F = F1*F2;
    
        
    case 5 % Gamma Distribution mu = fact*lambda
        fact = mu/lambda;
        % for lambda
        a = lambda^2/dlambda^2;
        invb = lambda/dlambda^2;
    
        % integral over gamma distribution
        fun1 = @(x) 1./(x+2*fact.*x).*1./gamma(a).*invb.^a.*x.^(a-1).*exp(-invb.*x);
        F1 = integral(fun1, 0, Inf);
        F1 = 1/F1;
        
        fun2 = @(x) x./(x+2*fact.*x).*1./gamma(a).*invb.^a.*x.^(a-1).*exp(-invb.*x);
        
        F2 = integral(fun2, 0, Inf);
        F = F1*F2;
        
    
        
     case 6  % Uniform Distribution mu = fact*lambda
         fact = mu/lambda;
        lower = (2*lambda-sqrt(12)*dlambda)/2; upper = (2*lambda+sqrt(12)*dlambda)/2; % for lambda
        F = (1/(1+2*fact)/(upper-lower)*log(upper/lower))^(-1)*1/(1+2*fact);
        
    case 7 % Uniform Distribution qratio involved
            fact = mu/lambda;
            [lambdaa, lambdab, mua, mub] = QratiotoLambda(qratio, epsl, fact, lambda);
            q = qratio*(mub - mua);
            if qratio == 0.5
                lambda1 = (lambdaa +lambdab)*0.5;
                lambda2 = lambda1;
            elseif qratio < 0.5
                lambda1 = (q+mua)/fact;
                lambda2 = (-q+mub)/fact;
            elseif qratio > 0.5
                lambda1 = (-q+mub)/fact;
                lambda2 = (q+mua)/fact;
            end
            fxy = -fact/q/(q+2*fact*(lambdaa-lambdab));    


            fun11 = @(y,x) 1./(y+2*x)*fxy;
            F11 = integral2(fun11, lambdaa, lambda1, mua, @(y)fact.*y+q );

            if qratio == 0.5
                F12 = 0;
            elseif qratio < 0.5
                fun12 = @(y,x) 1./(y+2*x)*fxy;
                F12 = integral2(fun12, lambda1, lambda2, @(y)fact.*y-q, @(y)fact.*y+q);
            elseif qratio > 0.5
                  fun12 = @(y,x) 1./(y+2*x)*fxy;
                F12 = integral2(fun12, lambda1, lambda2, mua, mub);
            end

            fun13 = @(y,x) 1./(y+2*x)*fxy;
            F13 = integral2(fun13, lambda2, lambdab, @(y)fact.*y-q, mub);

            F1 =  F11 + F12 + F13;
            F1 = 1/F1;


            fun21 = @(y,x) y./(y+2*x)*fxy;
            F21 = integral2(fun21, lambdaa, lambda1, mua, @(y)fact.*y+q );

            if qratio == 0.5
                F22 = 0;
            elseif qratio < 0.5
                fun22 = @(y,x) y./(y+2*x)*fxy;
                F22 = integral2(fun22, lambda1, lambda2, @(y)fact.*y-q, @(y)fact.*y+q);
            elseif qratio > 0.5
                fun22 = @(y,x) y./(y+2*x)*fxy;
                F22 = integral2(fun22, lambda1, lambda2, mua, mub);
            end

            fun23 = @(y,x) y./(y+2*x)*fxy;
            F23 = integral2(fun23, lambda2, lambdab, @(y)fact.*y-q, mub);

            F2 =  F21 + F22 + F23;
            F = F1*F2;

end