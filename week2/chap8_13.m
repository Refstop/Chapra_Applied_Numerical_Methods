% 8-13
clc; clear all; close all;
g = 9.81; k = 10;
m1 = 2; m2 = 3; m3 = 2.5;

A = [k 0; k -k];
b = [-m1*g; m2*g];
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
