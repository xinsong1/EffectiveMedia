function A = c11(lambda, mu, eta, epsl, epsm, alpha, factor,qratio)
% Closed form of effective medium component c11 = C1111 = A (in Backus notation)
% Assumes dlambda/lambda = dmu/mu = epsilon
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
            
% By Thomas Jordan 2011/03 for case 1
% modified by Xin Song 10/27/2014 11/06/2014 2/10/2016 2/21/2016 2/22/2016
% 2/23/2016 2/26/2016

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
          A = lambda + 2*mu...
                    + 4*(dlambda^2 + 2*dlm + dmu^2)*g11(lambda,mu,eta)...
                    + 4*(dlambda^2+dlm)*g13(lambda,mu,eta)...
                    + dlambda^2*g33(lambda,mu,eta)...
                    - 4*(dlambda^2+2*dlm)*g44(lambda,mu,eta);
    case 2 % Backus 2nd Order Approxiamte Theory
        A = (lambda + 2*mu)*...
            (1 - dlambda^2/(lambda+2*mu)^2);
    case 3 % Backus Exact Theory () lambda and mu are i.i.d
        % for muGamma
        ax = mu^2/dmu^2;
        invbx = mu/dmu^2;
        % for lambda
        ay = lambda^2/dlambda^2;
        invby = lambda/dlambda^2;
     
        % integral over gamma distribution
        % <a-f^2c^-1>
        fun1 = @(x,y) ((y+2.*x) - y.^2./(y+2.*x)).*1./gamma(ax).*invbx.^ax.*x.^(ax-1).*exp(-invbx.*x) ...
            .*1./gamma(ay).*invby.^ay.*y.^(ay-1).*exp(-invby.*y);
         A1 = integral2(fun1, 0, Inf, 0, Inf);
        
        % <c^-1>^-1
        fun2 = @(x,y) 1./(y+2*x).*1./gamma(ax).*invbx.^ax.*x.^(ax-1).*exp(-invbx.*x) ...
            .*1./gamma(ay).*invby.^ay.*y.^(ay-1).*exp(-invby.*y);
        A2 = integral2(fun2, 0, Inf, 0, Inf);
        A2 = 1/A2;
        
        % <fc^-1)^2
        fun3 = @(x,y) y./(y+2*x).*1./gamma(ax).*invbx.^ax.*x.^(ax-1).*exp(-invbx.*x) ...
            .*1./gamma(ay).*invby.^ay.*y.^(ay-1).*exp(-invby.*y);
        A3 = integral2(fun3, 0, Inf, 0, Inf);
        A3 = A3^2;
        A = A1+A2*A3;
        

    case 4 % Backus Exact Theory (Uniform) lambda and mu are i.i.d

        lowerx = (2*mu-sqrt(12)*dmu)/2; upperx = (2*mu+sqrt(12)*dmu)/2; %for mu
        lowery = (2*lambda-sqrt(12)*dlambda)/2; uppery = (2*lambda+sqrt(12)*dlambda)/2; % for lambda
        
        % integral over uniform distribution
        % <a-f^2c^-1>
        fun1 = @(x, y) ((y+2.*x) - y.^2./(y+2.*x)).*(upperx-lowerx)^-1.*(uppery-lowery)^-1;
        A1 = integral2(fun1, lowerx, upperx, lowery, uppery);
        
        % <c^-1>^-1
        fun2 = @(x,y) 1./(y+2*x).*(upperx-lowerx)^-1.*(uppery-lowery)^-1;
        A2 = integral2(fun2, lowerx, upperx, lowery, uppery);
        A2 = 1/A2;
        
        % <fc^-1)^2
        fun3 = @(x,y) y./(y+2*x).*(upperx-lowerx)^-1.*(uppery-lowery)^-1;
        A3 = integral2(fun3,lowerx, upperx, lowery, uppery);
        A3 = A3^2;
        A = A1+A2*A3;
        
    
    case 5 % Gamma Distribution mu = fact*lambda
        fact = mu/lambda;

        % for lambda
        a = lambda^2/dlambda^2;
        invb = lambda/dlambda^2;
     
        % integral over gamma distribution
        % <a-f^2c^-1>^-1
        fun1 = @(x) ((x+2*fact.*x) - x.^2./(x+2*fact.*x)).*1./gamma(a).*invb.^a.*x.^(a-1).*exp(-invb.*x);
        A1 = integral(fun1, 0, Inf);
        
        % <c^-1>^-1
        fun2 = @(x) 1./(x+2*fact.*x).*1./gamma(a).*invb.^a.*x.^(a-1).*exp(-invb.*x);
        A2 = integral(fun2, 0, Inf);
        A2 = 1/A2;
        
        % <fc^-1)^2
        fun3 = @(x) x./(x+2*fact.*x).*1./gamma(a).*invb.^a.*x.^(a-1).*exp(-invb.*x);
        A3 = integral(fun3, 0, Inf);
        A3 = A3^2;
        A = A1+A2*A3;    
        
        
     case 6 % Uniform Distribution mu = fact*lambda
       fact = mu/lambda;
       lower = (2*lambda-sqrt(12)*dlambda)/2; upper = (2*lambda+sqrt(12)*dlambda)/2; % for lambda
       
       A = (4*fact+4*fact^2)/(1+2*fact)/2*(upper+lower)+(1/(1+2*fact)/(upper-lower)*log(upper/lower))^(-1)*(1/(1+2*fact))^2;
       
       
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
        % integral over uniform distribution
        % <a-f^2c^-1>^-1
        fun11 = @(y, x) ((y+2.*x) - y.^2./(y+2.*x))*fxy;
        A11 = integral2(fun11,lambdaa, lambda1, mua, @(y)fact.*y+q );



        if qratio ==0.5

            A12 = 0;
        elseif qratio < 0.5
            fun12 = @(y, x) ((y+2.*x) - y.^2./(y+2.*x))*fxy;
            A12 = integral2(fun12, lambda1, lambda2, @(y)fact.*y-q, @(y)fact.*y+q);
        elseif qratio > 0.5
            fun12 = @(y, x) ((y+2.*x) - y.^2./(y+2.*x))*fxy;
            A12 = integral2(fun12, lambda1, lambda2, mua, mub);
        end 

        fun13 = @(y, x) ((y+2.*x) - y.^2./(y+2.*x))*fxy;
        A13 = integral2(fun13, lambda2, lambdab, @(y)fact.*y-q, mub);

        A1 = A11 +A12 + A13;

        % <c^-1>^-1
        fun21 = @(y, x) 1./(y+2*x)*fxy;
        A21 = integral2(fun21, lambdaa, lambda1, mua, @(y)fact.*y+q );

        if qratio == 0.5
            A22 = 0;
        elseif qratio < 0.5
            fun22 = @(y, x) 1./(y+2*x)*fxy;
            A22 = integral2(fun22, lambda1, lambda2, @(y)fact.*y-q, @(y)fact.*y+q);
        elseif qratio > 0.5
            fun22 = @(y, x) 1./(y+2*x)*fxy;
            A22 = integral2(fun22, lambda1, lambda2, mua, mub);

        end

        fun23 = @(y, x) 1./(y+2*x)*fxy;
        A23 = integral2(fun23, lambda2, lambdab, @(y)fact.*y-q, mub);

        A2 = A21 + A22 + A23;
        A2 = 1/A2;

        % <fc^-1)^2
        fun31 = @(y, x) y./(y+2*x)*fxy;
        A31 = integral2(fun31,lambdaa, lambda1, mua, @(y)fact.*y +q );

        if qratio == 0.5
            A32 = 0;
        elseif qratio < 0.5
            fun32 = @(y, x) y./(y+2*x)*fxy;
            A32 = integral2(fun32, lambda1, lambda2, @(y)fact.*y-q, @(y)fact.*y+q);

        elseif qratio > 0.5
            fun32 = @(y, x) y./(y+2*x)*fxy;
            A32 = integral2(fun32,lambda1, lambda2, mua, mub);

        end 
        fun33 = @(y, x) y./(y+2*x)*fxy;
        A33 = integral2(fun33, lambda2, lambdab, @(y)fact.*y-q, mub);

        A3 = A31 + A32 + A33;
        A3 = A3^2;

        A = A1+A2*A3;

   
end
end