% test function 
% plot 2 figures:
% (1) recaled kneer tensor v.s. eta
% (2) effective Cij v.s. eta
% Xin Song 01/10/2017

clear; close all; clc

plotpth = './'; % you may define your own plot path
%% Global parameters
% elastic parameters
lambda = 5.9e9;
mu = 5.9e9;

% isotropic case
c11o = lambda + 2*mu;
c12o = lambda;
c33o = lambda + 2*mu;
c13o = lambda;
c44o = mu;
c55o = mu;


%% Recaled Kneer tensor plot
x = logspace(-3.5,3.5, 50);
green11 = zeros(50,1);
green33 = zeros(50,1);
green13 = zeros(50,1);
green12 = zeros(50,1);
green44 = zeros(50,1);
green55 = zeros(50,1);
Tint = zeros(50, 1);
 ETA2 = zeros(50,1);
for n = 1:50
%     eta = 1 + x(n);
   
    eta = x(n);
    ETA2(n) = eta;
    green11(n) = g11(lambda,mu,eta);
    green33(n) = g33(lambda,mu,eta);
    green12(n) = g12(lambda,mu,eta);
    green13(n) = g13(lambda,mu,eta);
    green44(n) = g44(lambda,mu,eta);
    green55(n) = g55(lambda,mu,eta);
    
    Tint(n, 1) = tint2(eta^2-1) - tint4(eta^2-1); 
end

figure;
x = log10(x);
plot(x,green11,x,green33,x,green13, x, green12, x,green44,x,green55, 'Linewidth', 2)
legend('green11', 'green33', 'green13', 'green12', 'green44', 'green55', 'location', 'SouthWest')
set(gca, 'fontsize', 12)
xlabel('log(\eta)', 'fontsize', 18)
ylabel('g', 'fontsize', 18)
grid on
xlim([-3 3])
name = [plotpth, 'g-eta'];
print('-r400','-dpdf', name)





%% C v.s. eta
N = 500;
epsilon = .2;
x = logspace(-3,3, N);
ceta11 = zeros(N,1);
ceta33 = zeros(N,1);
ceta12 = zeros(N,1);
ceta13 = zeros(N,1);
ceta44 = zeros(N,1);
ceta55 = zeros(N,1);

Eta = sqrt(x);

for n = 1:N

    eta = sqrt(x(n));
    ceta11(n) = c11(lambda,mu,eta,epsilon,epsilon, 1, 1)/c11o-1;
    ceta33(n) = c33(lambda,mu,eta,epsilon,epsilon, 1, 1)/c33o-1;
    ceta12(n) = c12(lambda,mu,eta,epsilon,epsilon, 1, 1)/c12o-1;
    ceta13(n) = c13(lambda,mu,eta,epsilon,epsilon, 1, 1)/c13o-1;
    ceta44(n) = c44(lambda,mu,eta,epsilon,epsilon, 1, 1)/c44o-1;
    ceta55(n) = c55(lambda,mu,eta,epsilon,epsilon, 1, 1)/c55o-1;
end

x = log10(x);

figure;
plot(x,ceta11,x,ceta33,x,ceta13,x,ceta44,x,ceta55, 'linewidth', 2)
legend('C11', 'C33', 'C13', 'C44', 'C55', 'location', 'SouthWest')
set(gca, 'fontsize', 14)
hold on
plot([0,0], [0, -0.04],'k--', 'linewidth', 2)
xlabel('log(\eta^2)', 'fontsize', 18);
ylabel('\deltac/c', 'fontsize', 18)
title('\epsilon = 0.2', 'fontsize', 18)
text(-2,-0.0025, '\eta < 1', 'HorizontalAlignment', 'center', 'FontSize', 14)
text(1,-0.0025, '\eta > 1', 'HorizontalAlignment', 'center', 'FontSize', 14)
grid on
name = [plotpth, 'c-eta' ];
% print('-r200', '-djpeg', name)




