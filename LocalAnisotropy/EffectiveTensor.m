function [An, Cn, Bn, Fn, Nn, Ln, A, C, B, F, N, L] = EffectiveTensor(C, xi, eta)
% this function calculated 6 TI component Effective Tensor for locally
% anisotropy
% Input
% local stiffness tensor with maximum hexgonal symmetry (cubic, isotropic work) 
% C   = [[C1111  C1122 C1133      0      0      0]; ...
%            [C1122  C1111  C1133      0      0      0]; ...
%            [C1133   C1133 C3333      0      0      0]; ...
%            [    0      0      0   C1313      0      0]; ...
%            [    0      0      0       0  C1313      0]; ...
%            [    0      0      0       0      0  C1212]]; 
% xi: 0 < xi < inf 
% eta: 0 < eta < inf

% Example
% C   = [[217.22  71.66  68.54      0      0      0]; ...
%            [71.66  217.22  68.54      0      0      0]; ...
%            [68.54   68.54 243.47      0      0      0]; ...
%            [    0      0      0   79.32      0      0]; ...
%            [    0      0      0       0  79.32      0]; ...
%            [    0      0      0       0      0  72.78]]; 
%        
%        
% xi = 1;
% eta = 0.1;


% Output
% An, Cn, Bn, Fn, Ln are the effective values with second order
% perturbation
% A, C, B, F, N, L are the Voigt average

% Xin Song 2016/06

Ai = C(1,1);
Ci = C(3,3);
Fi = C(1,3);
Ni = C(6,6);
Li = C(4,4);
Bi = Ai-2*Ni;

a = Ai - 2*Ni;
b = Ni;
r = Ai + Ci - 2*Fi - 4*Li;
d = -Ai + Fi+ 2*Ni;
e = Li - Ni;

A = mA(xi, a, b, r, d, e);
B = mB(xi, a, r, d);
F = mF(xi, a, r, d);
C = mC(xi, a, b, r, d, e);
N = mN(xi, b, r, e);
L = mL(xi, b, r, e);

AA = mAA(xi, a, b, r, d, e);
BB = mBB(xi, a, r, d);
FF = mFF(xi, a, r, d);
CC = mCC(xi, a, b, r, d, e);
NN = mNN(xi, b, r, e);
LL = mLL(xi, b, r, e);
AB = mAB(xi, a, b, r, d, e);
AF = mAF(xi, a, b, r, d, e);
AC = mAC(xi, a, b, r, d, e);
BF = mBF(xi, a, r, d);
FD = mFD(xi, a, r, d);
CB = mCB(xi, a, b, r, d, e);
CF = mCF(xi, a, b, r, d, e);

AH = mAH(xi, a, b, r, d, e);
AD = mAD(xi, a, b, r, d, e);
BD = mBD(xi, a, r, d);
HF = mHF(xi, a, b, r, d, e);

C14 = msC14(xi, r, d, e);
C15 = msC15(xi, r, d, e);
C16 = msC16(xi, r, d);

C34 = msC34(xi, r, d);
C35 = msC35(xi, r, d, e);
C36 = msC36(xi, r, d, e);

C45 = msC45(xi, r, e);
C46 = msC46(xi, r, e);
C56 = msC56(xi, r, e);

C1442 = mC1442(xi, r, d, e);
C1552 = mC1552(xi, r, d, e);
C1662 = mC1662(xi, r, d, e);

C1443 = mC1443(xi, r, d, e);
C1553 = mC1553(xi, r, d, e);
C1663 = mC1663(xi, r, d, e);


% calculate the Kneer tensor

E11 = KneerE11(eta, A, C, F, L, N);
E12 = KneerE12(eta, A, C, F, L, N);
E13 = KneerE13(eta, A, C, F, L);
E33 = KneerE33(eta, A, C, F, L);
E44 = KneerE44(eta, A, C, F, L, N);
E55 = KneerE55(eta, A, C, F, L, N);



dA = E11.*(AA - A.*A) + E11.*(BB - B.*B) ...
    + 2*E12.*(AB - A.*B) ...
    + 2*E13.*(AF - A.*F) ...
    + 2*E13.*(BF - B.*F) ...
    + E33.*(FF - F.*F) ...
    + 4*E44.*C14 + 4*E55.*C15 + 4*E55.*C16;

dC = 2* E11.*(FF - F.*F) ...
    + 2*E12.*(FD - F.*F) ...
    + 4*E13.*(CF - C.*F) ...
    + E33.*(CC - C.*C) ...
    + 4*E44.*C34 + 4*E55.*C35 + 4*E55.*C36;

dB = 2* E11.*(AB - A.*B) ...
    + E12.*(AH - A.*A)  + E12.*(BB - B.*B)...
    + E13.*(AD - A.*F) + E13.*(BF - B.*F)...
    + E13.*(BD - B.*F) + E13.*(HF - A.*F) ...
    + E33.*(FD - F.*F) ...
    + 4*E44.*C1442 + 4*E55.*C1552 + 4*E55.*C1662;

dF = E11.*(AF - A.*F) + E11.*(BF - B.*F) ...
    + E12.*(AD - A.*F) + E12.*(BF - B.*F) ...
    + E13.*(FF - F.*F) + E13.*(FD - F.*F)...
    + E13.*(AC - A.*C) + E13.*(CB - C.*B) ...
    + E33.*(CF - C.*F) ...
    + 4*E44.*C1443 + 4*E55.*C1553 + 4*E55.*C1663;


dN = 4*E44.*(NN - N.*N) ...
    + 2*E11.*C14 + 2*E12.*C1442 ...
    + 4*E13.*C1443 + E33.*C34 ...
    + 4*E55.*C45 + 4*E55.*C46;


dL = 4*E55.*(LL - L.*L) ...
    + E11.*C15 + E11.*C16 + 2*E12.*C1552 ...
    + 2*E13.*C1553  + 2*E13.*C1663 ...
    + E33.*C35 ...
    + 4*E44.*C45 + 4*E55.*C56;
% new elastic tensor C bar + deltaC
An = dA + A;
Cn = dC + C;
Bn = dB + B;
Fn = dF + F;
Nn = dN + N;
Ln = dL + L;
