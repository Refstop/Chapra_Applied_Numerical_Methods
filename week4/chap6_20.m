clc; clear all; close all;
k1 = 40000; k2 = 40; m = 95; g = 9.81; h = 0.43;

%% 연속대입법
x = 0; % 초기값
es = 0.01;
F = @(d) (2*k2*d.^2.5)/(5*m*g) + (k1*d.^2)/(2*m*g) - h;
i = 0;
while(1)
    xr = F(x);
    ea = abs((xr-x)/xr)*100;
    x = xr;
    if ea <= es | i>=100
        break;
    end
    i=i+1;
end
if i==100
    fprintf("발산하여 근을 구할 수 없음.\n");
else
    fprintf("연속대입법: %.10f\n", xr);
end

%% Newton-Raphson
es=0.0001;
i=0;
xr=5;
F = @(d) (2*k2*d^2.5)/5 + 0.5*k1*d^2 - m*g*d - m*g*h;
dF = @(d) k2*d^1.5 + k1*d - m*g;
while(1)
    xrold=xr;
    xr = xr-F(xr)/dF(xr);
    i=i+1;
    if xr~=0
        ea=abs((xr-xrold)/xr)*100;
    end
    if ea<=es | i>=100
        break;
    end
end
fprintf("Newton Raphson: %.10f\n", xr);


%% 할선법
es=0.01;
i=0;
xrold=1; xr=5;
F = @(d) (2*k2*d^2.5)/5 + 0.5*k1*d^2 - m*g*d - m*g*h;
while(1)
    xroldold = xrold;
    xrold = xr;
    xr = xrold-(F(xrold)*(xroldold-xrold))/(F(xroldold)-F(xrold));
    i=i+1;
    if xr~=0
        ea=abs((xr-xrold)/xr)*100;
    end
    if ea<=es | i>=100
        break;
    end
end
fprintf("할선법: %.10f\n", xr);

%% 수정 할선법
es=0.01;
i=0;
xr=5;
delta_xr=0.01;
F = @(d) (2*k2*d^2.5)/5 + 0.5*k1*d^2 - m*g*d - m*g*h;
while(1)
    xrold = xr;
    xr = xr-(delta_xr*F(xr))/(F(xr+delta_xr)-F(xr));
    i=i+1;
    if xr~=0
        ea=abs((xr-xrold)/xr)*100;
    end
    if ea<=es | i>=100
        break;
    end
end
fprintf("수정 할선법: %.10f\n", xr);

%% Brent법
F = @(d) (2*k2*d^2.5)/5 + 0.5*k1*d^2 - m*g*d - m*g*h;
a = 0; b = 5;
x1 = a; x2 = b; x3 = (x1+x2)/2;
f1 = F(x1); f2 = F(x2);
tol = 0.01;
for k=1:30
    f3=F(x3);
    if abs(f3)<tol
        xr=x3;
        break;
    end
    if f1*f3<0.0; b=x3; % 근이 추정근보다 작으면 (a, b=x3, b)
    else, a=x3; % 근이 추정근보다 크면 (a, a=x3, b)
    end
    if (b-a)<tol*max(abs(b),1.0)
        xr=0.5*(a+b);
        break;
    end
    % Second secant method
    denom=(f2-f1)*(f3-f1)*(f2-f3);
    numer=x3*(f1-f2)*(f2-f3+f1)+f2*x1*(f2-f3)+f1*x2*(f3-f1);
    if denom==0, dx=b-a;
    else dx=f3*numer/denom;
    end
    x=x3+dx; % 역 2차 보간법에 따른 추정근
    if (b-x)*(x-a)<0.0 % 추정근이 범위를 벗어나면 이분법으로
        dx=0.5*(b-a); x=a+dx; % 이분법에 따른 새로운 추정근
    end
    if x<x3 % 추정근이 전 추정값보다 작으면
        x2=x3; f2;f3; % (x1, x2=x3, x2)
    else % 크면
        x1=x3; f1=f3; % (x1, x1=x3, x2)
    end
    x3=x; % (x1, x3=x, x2)
end
fprintf("Brent법: %.10f\n", xr);