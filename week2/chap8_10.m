% 8-10
clc; clear all; close all;
% F1 F2 F3 H2 V2 V3
% Gauss
A = [cosd(30) 0 -cosd(60) 0 0 0;
    0 1 cosd(60) 0 0 0;
    sind(30) 0 sind(60) 0 0 0;
    cosd(30) 1 0 1 0 0;
    sind(30) 0 0 0 1 0;
    0 0 sind(60) 0 0 1];
b = [1000; 0; 1000; 0; 0; 0];
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