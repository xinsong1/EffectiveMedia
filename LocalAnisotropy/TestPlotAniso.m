% test function 
% plot 2 subplots:
% (1) voigt average v.s. xi
% (2) effective Cij v.s. xi
% Xin Song 01/11/2017
clear; clc; close all


plotpth = './'; % define your own plotpth


%% tensor 

C_oliv = [[212.5  79.3  69.9     0     0     0]; ...
            [79.3 212.4  69.9     0     0     0]; ...
            [69.9  69.9 320.5     0     0     0]; ...
            [0     0     0  77.9     0     0]; ...
            [0     0     0     0  77.9     0]; ...
            [0     0     0     0     0  66.5]];
        

name = 'Oliv';

Tensor = C_oliv;


Ai = Tensor(1,1);
Ci = Tensor(3,3);
Fi = Tensor(1,3);
Ni = Tensor(6,6);
Li = Tensor(4,4);
Bi = Ai-2*Ni;

a = Ai - 2*Ni;
b = Ni;
r = Ai + Ci - 2*Fi - 4*Li;
d = -Ai + Fi+ 2*Ni;
e = Li - Ni;

%% xi from 10^-3 10^-3
Xi = logspace(-3,3);
Eta = [0.03, 1, 30];
n = length(Xi);
m = length(Eta);

A = zeros(n, m);
B = zeros(n, m);
F = zeros(n, m);
C = zeros(n, m);
N = zeros(n, m);
L = zeros(n, m);


An = zeros(n, m);
Bn = zeros(n, m);
Fn = zeros(n, m);
Cn = zeros(n, m);
Nn = zeros(n, m);
Ln = zeros(n, m);

for j = 1:m
    eta = Eta(j);
    for i = 1:n

        [Ann, Cnn, Bnn, Fnn, Nnn, Lnn, AA, CC, BB, FF, NN, LL] = EffectiveTensor(Tensor, Xi(i), eta);

        % Voigt Average
        A(i, j) = AA;
        C(i, j) = CC;
        B(i, j) = BB;
        F(i, j) = FF;
        N(i, j) = NN;
        L(i, j) = LL;

        An(i, j) = Ann;
        Cn(i, j) = Cnn;
        Bn(i, j) = Bnn;
        Fn(i, j) = Fnn;
        Nn(i, j) = Nnn;
        Ln(i ,j) = Lnn;

    end
end

% values when xi =1
ai = mA(1, a, b, r, d, e);
ci = mC(1, a, b, r, d, e);
ni = mN(1, b, r, e);
li = mL(1, b, r, e);
fi = mF(1, a, r, d);
bi = mB(1, a, r, d);

Xi = log10(Xi);


colmap=[255, 51, 51; 255, 153, 51; 0, 102, 0;0, 204, 204;  51, 51, 255 ; 153, 51, 255];
colmap = colmap/255;


%% figure11 for the paper
figure;
pos = [0.05, 0.3, 0.4, 0.4];
subplot('position', pos)
l= 1; % the voigt average is same for different l bc voigt average is not varying with xsi
plot(Xi, A(:, l), 'Color', colmap(1, :), 'linewidth', 2) 
hold on
plot(Xi, B(:, l), 'Color', colmap(2, :), 'linewidth', 2)
hold on
plot(Xi, F(:, l), 'Color', colmap(3, :), 'linewidth', 2) 
hold on
plot(Xi, C(:, l), 'Color', colmap(4, :), 'linewidth', 2)
hold on
plot(Xi, N(:, l), 'Color', colmap(5, :), 'linewidth', 2) 
hold on
plot(Xi, L(:, l), 'Color', colmap(6, :), 'linewidth', 2)
hold on
plot(Xi(1), Ai, 'ko', Xi(1), Bi, 'ko', Xi(1), Fi, 'ko', Xi(1), Ci, 'ko', ...
    Xi(1), Ni, 'ko', Xi(1), Li, 'ko', 'Markersize', 5, 'MarkerFaceColor', [0, 0, 0])
 xlim([-3, 3])
set(gca, 'fontsize', 10)
grid on
xlabel('log_{10}(\xi)', 'fontsize', 10)
title('Voigt Average')

pos = [0.55, 0.3, 0.4, 0.4];
subplot('position', pos)
l= 1; % the voigt average is same for different l bc voigt average is not varying with xsi
plot(Xi, (A(:, l)-ai)/ai, 'Color', colmap(1, :), 'linewidth', 2) 
hold on
plot(Xi, (B(:, l)-bi)/bi, 'Color', colmap(2, :), 'linewidth', 2)
hold on
plot(Xi, (F(:, l)-fi)/fi, 'Color', colmap(3, :), 'linewidth', 2) 
hold on
plot(Xi, (C(:, l)-ci)/ci, 'Color', colmap(4, :), 'linewidth', 2)
hold on
plot(Xi, (N(:, l)-ni)/ni, 'Color', colmap(5, :), 'linewidth', 2) 
hold on
plot(Xi, (L(:, l)-li)/li, 'Color', colmap(6, :), 'linewidth', 2)
hold on
legend('C_{11}', 'C_{12}', 'C_{13}', 'C_{33}', 'C_{44}', 'C_{55}')
plot(Xi(1), (Ai-ai)/ai, 'ko', Xi(1), (Bi-bi)/bi, 'ko', Xi(1), (Fi-fi)/fi, 'ko', ...
    Xi(1), (Ci-ci)/ci, 'ko', Xi(1), (Ni-ni)/ni, 'ko', Xi(1), (Li-li)/li, 'ko', ...
    'Markersize', 5, 'MarkerFaceColor', [0, 0, 0])
 xlim([-3, 3])
set(gca, 'fontsize', 10)
grid on
xlabel('log_{10}(\xi)', 'fontsize', 10)
title('Effective C^*')
print('-r200', '-dpdf', [plotpth, 'C-xi' ])


