clc;clear all;close all;
Fx1 = @(x2, x3) (20+x2-2*x3)/8;
Fx2 = @(x1, x3) -(-38-2*x1+x3)/6;
Fx3 = @(x1, x2) (-34+3*x1+x2)/7;
lambda=1.2;

%% 역함수
A = [2 -6 -1; -3 -1 7; -8 1 -2];
b = [-38; -34; -20];
X = inv(A)*b;
fprintf("역함수: %.8f %.8f %.8f\n", X(1), X(2), X(3));


%% 가우스 소거법
A = [2 -6 -1; -3 -1 7; -8 1 -2];
b = [-38; -34; -20];
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
fprintf("가우스 소거법(Gauss Elimination): %.8f %.8f %.8f\n", x(1), x(2), x(3));

%% 연속대입법(Gauss Seidal)
xr1 = 0; xr2 = 0; xr3 = 0;
es = 0.001; i=0;
while(i<=1000)
    xr1old = xr1;
    xr1 = Fx1(xr2, xr3);
    xr2 = Fx2(xr1, xr3);
    xr3 = Fx3(xr1, xr2);
    ea = abs((xr1-xr1old)/xr1)*100;
    if ea <= es
        break;
    end
    i=i+1;
end
if i>=1000
    fprintf("방정식의 기울기 값의 절댓값이 1보다 큰 경우가 존재하므로 발산한다.\n");
else
    fprintf("연속대입법(Gauss Seidal): %.8f %.8f %.8f\n", xr1, xr2, xr3);
end


%% 연속대입법(Jacobi)
xr1 = 0; xr2 = 0; xr3 = 0;
es = 0.001; i=0;
while(i<=1000)
    xr1old = xr1;
    xr2old = xr2;
    xr3old = xr3;
    xr1 = Fx1(xr2old, xr3old);
    xr2 = Fx2(xr1old, xr3old);
    xr3 = Fx3(xr1old, xr2old);
    ea = abs((xr1-xr1old)/xr1)*100;
    if ea <= es
        break;
    end
    i=i+1;
end
if i>=1000
    fprintf("방정식의 기울기 값의 절댓값이 1보다 큰 경우가 존재하므로 발산한다.\n");
else
    fprintf("연속대입법(Jacobi): %.8f %.8f %.8f\n", xr1, xr2, xr3);
end
