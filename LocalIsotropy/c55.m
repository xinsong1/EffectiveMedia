function L = c55(lambda, mu, eta, epsl, epsm, alpha,  factor, pratio)
% Closed form of effective medium component c55 = C1313
% Assumes dlambda/lambda = dmu/mu = epsilon
% Elastic parameter calculation:
            % Factor = 1----------Jordan's EM 
            % Factor = 2----------Backus Approximate Theory
            % Factor = 3----------Backus Exact Theory (Gamma)
            % Factor = 4----------Backus Exact Theory (Uniform)
            % Factor = 5----------Backus Exact Theory (Gamma lambda = fact*mu)
            % Factor = 6----------Backus Exact Theory (Uniform lambda = fact*mu)
            % Factor = 7----------Backus exact Theory (Uniform, qratio involved)
% By Thomas Jordan 2011/03 for case 1
% modified by Xin Song 10/27/2014 11/06/2014 2/21/2016 2/22/2016

if lambda <= 0, error('lambda <= 1'), end
if mu <= 0, error('mu <= 1'), end
if epsl < 0, error('epsl <= 0'), end
if epsm < 0, error('epsm <= 0'), end
% if eta <= 1, error('eta <= 1'), end

dmu = epsm*mu;

switch factor
    case 1 % Jordan Effective Medium Theory
            L = mu...
                + 4*dmu^2*g55(lambda,mu,eta);
    case 2 % Backus Approximate Theory
            L = mu*...
                (1 - dmu^2/mu^2);
    case 3 % Backus Exact Theory (Gamma, lambda and mu are iid)
            a = mu^2/dmu^2;
            invb = mu/dmu^2; 
    
            % integral over gamma distribution
            fun = @(x) 1./x.*1./gamma(a).*invb.^a.*x.^(a-1).*exp(-invb.*x);
            L = integral(fun, 0, Inf);
            L = 1/L;
            
    case 4 % Backus Exact Theory (Uniform, lambda and mu are iid)
        lower = (2*mu-sqrt(12)*dmu)/2; upper = (2*mu+sqrt(12)*dmu)/2;
        fun_U = @(x) 1./x.*(upper-lower)^-1;
        L = integral(fun_U, lower, upper);
        L = 1/L;
        
   
        
    case 5 % Backus Exact Theory (Gamma,  mu = fact*lambda)
        a = mu^2/dmu^2;
            invb = mu/dmu^2; 
    
            % integral over gamma distribution
            fun = @(x) 1./x.*1./gamma(a).*invb.^a.*x.^(a-1).*exp(-invb.*x);
            L = integral(fun, 0, Inf);
            L = 1/L;

     case 6 % Backus Exact Theory (Uniform,  mu = fact*lambda)
        lower = (2*mu-sqrt(12)*dmu)/2; upper = (2*mu+sqrt(12)*dmu)/2;
        
        L = 1/(upper-lower)*log(upper/lower);
        L = 1/L;
        
    case 7 % Uniform Distribution qratio involved
        
        lower = (2*mu-sqrt(12)*dmu)/2; upper = (2*mu+sqrt(12)*dmu)/2;
        fun_U = @(x) 1./x.*(upper-lower)^-1;
        L = integral(fun_U, lower, upper);
        L = 1/L;
  
end

end