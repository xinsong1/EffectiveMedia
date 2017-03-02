function C = c33(lambda, mu, eta, epsl, epsm, alpha, factor, qratio)
% Closed form of effective medium component c33 = C3333
% Assumes dlambda/lambda = dmu/mu% Elastic parameter calculation:
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
             
% By Thomas Jordan 2011/03
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
%         C = lambda + 2*mu...
%             + 4*dlambda^2*g11(lambda,mu,eta)...
%             + 4*dlambda*(dlambda + 2*dmu)*g13(lambda,mu,eta)...
%             + (dlambda + 2*dmu)^2*g33(lambda,mu,eta)...
%             - 4*dlambda^2*g44(lambda,mu,eta);
        
           C = lambda + 2*mu...
            + 4*dlambda^2*g11(lambda,mu,eta)...
            + 4*(dlambda^2 + 2*dlm)*g13(lambda,mu,eta)...
            + (dlambda^2 + 4*dmu^2 + 4*dlm)*g33(lambda,mu,eta)...
            - 4*dlambda^2*g44(lambda,mu,eta);
    case 2  % Backus 2nd Order Approxiamte Theory
%         C = (lambda + 2*mu)*...
%             (1 - (dlambda^2+4*dmu^2 + 4*dlm)/(lambda+2*mu)^2);
          C = lambda + 2*mu - (dlambda^2 + 4*dmu^2 + 4*dlm)/(lambda+2*mu);
    case 3 % Backus Exact Theory (Gamma) lambda and mu are i.i.d
        % for mu
        ax = mu^2/dmu^2;
        invbx = mu/dmu^2;
        % for lambda
        ay = lambda^2/dlambda^2;
        invby = lambda/dlambda^2;
    
        % integral over gamma distribution
        fun = @(x,y) 1./(y+2*x).*1./gamma(ax).*invbx.^ax.*x.^(ax-1).*exp(-invbx.*x) ...
            .*1./gamma(ay).*invby.^ay.*y.^(ay-1).*exp(-invby.*y);
        C = integral2(fun, 0, Inf, 0, Inf);
        C = 1/C;
    case 4 % Backus Exact Theory (Uniform) lambda and mu are i.i.d
        lowerx = (2*mu-sqrt(12)*dmu)/2; upperx = (2*mu+sqrt(12)*dmu)/2; %for mu
        lowery = (2*lambda-sqrt(12)*dlambda)/2; uppery = (2*lambda+sqrt(12)*dlambda)/2; % for lambda
        
        % integral over gamma distribution
        fun = @(x,y) 1./(y+2*x).*(upperx-lowerx)^-1.*(uppery-lowery)^-1;
        C = integral2(fun, lowerx, upperx, lowery, uppery);
        C = 1/C;
    
        
    case 5  % Gamma Distribution mu = fact*lambda
        fact = mu/lambda;
        % for lambda
        a = lambda^2/dlambda^2;
        invb = lambda/dlambda^2;
    
        % integral over gamma distribution
        fun = @(x) 1./(x+2*fact.*x).*1./gamma(a).*invb.^a.*x.^(a-1).*exp(-invb.*x);
        C = integral(fun, 0, Inf);
        C = 1/C;
        
    case 6 % Uniform Distribution mu = fact*lambda
          fact = mu/lambda;
        lower = (2*lambda-sqrt(12)*dlambda)/2; upper = (2*lambda+sqrt(12)*dlambda)/2; % for lambda
        
       C = 1/(1+2*fact)/(upper - lower)*log(upper/lower);
       C = 1/C;
        
       
    case 7   % Uniform Distribution qratio involved
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

        fun1 = @(y, x) 1./(y+2*x)*fxy;
        C1 = integral2(fun1, lambdaa, lambda1, mua, @(y)fact.*y +q);


        if qratio == 0.5
            C2 = 0;
        elseif qratio < 0.5
            fun2 = @(y, x) 1./(y+2*x)*fxy;
            C2 = integral2(fun2, lambda1, lambda2, @(y)fact.*y-q, @(y)fact.*y+q);
        elseif qratio > 0.5
            fun2 = @(y, x) 1./(y+2*x)*fxy;
            C2 = integral2(fun2, lambda1, lambda2, mua, mub);
        end

        fun3 = @(y, x) 1./(y+2*x)*fxy;
        C3 = integral2(fun3, lambda2, lambdab, @(y)fact.*y-q, mub);

        C = C1 + C2 + C3;
        C = 1/C;

  
end
end