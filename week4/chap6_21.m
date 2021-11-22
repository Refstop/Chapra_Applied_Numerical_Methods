clc;clear all;close all;
v0=30; x=90; y0=1.8; y=1; g=9.81;

%% 연속대입법
xr = 0; % 초기값
i = 0; es = 0.01; ea = 1;
F = @(theta) atan(g*x/(2*v0^2*cos(theta)^2) + (y-y0)/x);

while(ea>es)
    xrold = xr;
    xr = F(xr);
    if xr ~= 0
        ea = abs((xr-xrold)/xr)*100;
    end
    if i>=100
        break;
    end
    i=i+1;
end
fprintf("연속대입법: %.10f\n", rad2deg(xr));

%% Newton-Raphson
xr = 0.01;
i = 0; es = 0.0001; ea = 1;


F = @(th) tan(th)*x - g*x^2/(2*v0^2*cos(th)^2) + y0-y;
dF = @(th) x/tan(th) - g*x^2*tan(th)/(v0^2*cos(th)^2);
% F = @(d) exp(-d)-d;
% dF = @(d) -exp(-d)-1;
while(ea>es)
    xrold = xr;
    xr = xr-F(xr)/dF(xr);
    
    if xr~=0
        ea=abs((xr-xrold)/xr)*100;
    end
    if i>=100
        break;
    end
    i=i+1;
end
fprintf("Newton Raphson: %.10f\n", rad2deg(xr));


%% 할선법
i=2;
xr(i-1) = 0;
xr(i) = 0.01;
es = 0.01; ea = 1;

F = @(th) tan(th)*x - g*x^2/(2*v0^2*cos(th)^2) + y0-y;
while(ea>es)
    a = F(xr(i))*(xr(i-1)-xr(i));
    b = F(xr(i-1)) - F(xr(i));
    xr(i+1) = xr(i)-a/b;
    if xr(i+1)~=0
        ea=abs((xr(i+1)-xr(i))/xr(i+1))*100;
    end
    if i>=100
        break;
    end
    i=i+1;
end
fprintf("할선법: %.10f\n", rad2deg(xr(i)));

%% 수정 할선법
es=0.01; ea=1;
i=0;
xr=1;
delta=0.01;
F = @(th) tan(th)*x - g*x^2/(2*v0^2*cos(th)^2) + y0-y;
while(ea>es)
    xrold = xr;
    xr = xr-(delta*xr*F(xr))/(F(xr+delta*xr)-F(xr));
    i=i+1;
    if xr~=0
        ea=abs((xr-xrold)/xr)*100;
    end
    if i>=100
        break;
    end
end
fprintf("수정 할선법: %.10f\n", rad2deg(xr));

%% Brent법
F = @(th) tan(th)*x - g*x^2/(2*v0^2*cos(th)^2) + y0-y;
a = 0; b = 0.8;
x1 = a; x2 = b; x3 = (x1+x2)/2;
f1 = F(x1); f2 = F(x2);
tol = 0.01;
for k=1:30
    f3=F(x3);
    if abs(f3) < tol
        xr=x3;
        break;
    end
    if f1*f3 < 0.0
        b=x3; % 근이 추정근보다 작으면 (a, b=x3)
    else
        a=x3; % 근이 추정근보다 크면 (a=x3, b)
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
        x2=x3; f2=f3; % (x1, x2=x3, x2)
    else % 크면
        x1=x3; f1=f3; % (x1, x1=x3, x2)
    end
    x3=x; % (x1, x3=x, x2)
end
fprintf("Brent법: %.10f\n", rad2deg(xr));
