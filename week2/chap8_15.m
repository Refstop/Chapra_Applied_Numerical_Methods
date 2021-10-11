% 8-15
clc; clear all; close all;
R12 = 35; R23 = 30; R34 = 8; R45 = 15; R56 = 5; R25 = 10; R35 = 7;
V1 = 20; V6 = 140;
% i12, i23, i34, i45, i56, i25, i35
% all clockwise
A = [-R12 0 0 -R56 -R25 0;
    0 R23 0 0 R25 -R35;
    0 0 R34+R45 0 0 R35;
    0 0 1 1 -1 -1;
    1 1 0 0 -1 0;
    0 1 -1 0 0 1];
b = [-(V1-V6); 0; 0; 0; 0; 0];
[m, n] = size(A);
C = [A b];
for i = 1:m
    for j = (i+1):m
        C(j,:) = C(j,:) - (C(j,i)/C(i,i))*C(i,:);
    end
end

X_gauss = zeros(1,n)';

for i = m:-1:1
    X_gauss(i,1) = (C(i,n+1) - C(i,1:n)*X_gauss)/C(i,i);
end

X_gauss
% LU
[L,U] = lu(A);
Y = zeros(1,n)';
C_1 = [L b];
for i = 1:m
    Y(i,1) = (C_1(i,n+1) - C_1(i,1:n)*Y)/C_1(i,i);
end
C_2 = [U Y];
X_LU = zeros(1,n)';
for i = m:-1:1
    X_LU(i,1) = (C_2(i,n+1) - C_2(i,1:n)*X_LU)/C_2(i,i);
end
X_LU

% inverse
X_inv = inv(A)*b
