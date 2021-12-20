% 8-5
clc; clear all; close all;
%% 가우스 소거법
A = [3 + 2i 4; -i 1];
b = [2+i; 3];
[m, n] = size(A);
if m~=n
    fprintf("정방행렬이 아닙니다.\n");
end
% 전방소거법
nb = n+1;
Aug = [A b];
for k = 1:n-1
    for i = k+1:n
        factor = Aug(i,k)/Aug(k,k);
        Aug(i,k:nb) = Aug(i,k:nb) - factor*Aug(k,k:nb);
    end
end

% 후방대입법
x = zeros(n,1);
x(n) = Aug(n,nb)/Aug(n,n);
for i = n-1:-1:1
    x(i) = (Aug(i,nb)-Aug(i,i+1:n)*x(i+1:n))/Aug(i,i);
end
fprintf("가우스 소거법(Gauss Elimination): %.8f %.8f\n", x(1), x(2));

%% LU 분해법
% [L,U] = lu(A);
% 상삼각, 하삼각행렬 만들기
U = A;
L = zeros(n,n);
for i = 1:n-1
    L(i:n,i) = 1/U(i,i)*U(i:n,i);
    for j = i+1:n
        U(j,i:n) = U(j,i:n) - U(j,i)/U(i,i)*U(i,i:n);
    end
end
L(n,n) = 1/U(n,n)*U(n,n);

% 후방대입법
nb=n+1;
Y = zeros(n,1);
C_1 = [L b];
Y(1) = C_1(1,nb)/C_1(1,1);
for i = 2:n
    Y(i) = (C_1(i,nb)-C_1(i,1:i-1)*Y(1:i-1))/C_1(i,i);
end
X_LU = zeros(n,1);
C_2 = [U Y];
X_LU(n) = C_2(n,nb)/C_2(n,n);
for i = n-1:-1:1
    X_LU(i) = (C_2(i,nb)-C_2(i,i+1:n)*X_LU(i+1:n))/C_2(i,i);
end
fprintf("LU 분해법(LU Decomposition): %.8f %.8f\n", X_LU(1), X_LU(2));

%% inverse
X_inv = inv(A)*b;
fprintf("역행렬(Inverse Matrix): %.8f %.8f\n", X_inv(1), X_inv(2));