% 8-14
clc; clear all; close all;
R12 = 10; R23 = 20; R34 = 2; R45 = 5; R56 = 25; R25 = 5;
V1 = 150; V6 = 0;
% i12, i23(=i34=i45), i25, i56
% clockwise, counter-clockwise
A = [R12 0 R25 R56;
    0 R23+R34+R45 R25 0;
    1 1 -1 0;
    0 1 -1 1];
b = [-(V6-V1); 0; 0; 0];
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