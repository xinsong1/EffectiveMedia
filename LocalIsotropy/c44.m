function M = c44(lambda, mu, eta, epsl, epsm, alpha,  factor, pratio)
% Closed form of effective medium component c44 = C1212
% Assumes dlambda/lambda = dmu/mu = epsilon
% Elastic parameter calculation:
            % Factor = 1----------Jordan's EM 
            % Factor = 2----------Backus Approximate Theory
            % Factor = 3----------Backus Exact Theory (Gamma iid) 
            % Factor = 4----------Backus Exact Theory (Uniform iid)
            % Factor = 5----------Backus Exact Theory (Gamma mu = fact*lambda)
            % Factor = 6----------Backus Exact Theory (Uniform mu = fact*lambda)
            % Factor = 7----------Backus Exact Theory (Uniform, qratio involved)
% By Thomas Jordan 2011/03 for case 1
% modified by Xin Song 10/27/2014 2/21/2016 2/22/2016

if lambda <= 0, error('lambda <= 1'), end
if mu <= 0, error('mu <= 1'), end
if epsl < 0, error('epsl <= 0'), end
if epsm < 0, error('epsm <= 0'), end
% if eta <= 1, error('eta <= 1'), end
% dlambda = epsilon*lambda;

dmu = epsm*mu;
switch factor
    case 1
        M = mu...
            + 4*dmu^2*g44(lambda,mu,eta);
    case 2
        M = mu;
    case 3
        M = mu;
    case 4
        M = mu;
    case 5;
        M = mu;
    case 6
        M = mu;      
    case 7
        M = mu;
    
end
end