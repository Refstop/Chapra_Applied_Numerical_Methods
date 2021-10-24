clc;clear all;close all;
%% 역 2차 보간법
L=45; Ac=6.362*10^-4; E=1.5*10^11; W=9000;
F = @(d) ((2*Ac*E*d)/L)*(1 - L/sqrt(L^2+d^2)) - W;

x=0;

i=3;
Es=0.00001;

x(i-2)=0;
x(i-1)=3;
x(i)=4;


f(i-2) = F(x(i-2));
f(i-1) = F(x(i-1));
f(i) = F(x(i));

Ea(i-2) = 10;
Ea(i-1) = 10;
Ea(i) = 10;

while Ea(i) >= Es
    a = f(i-1)*f(i);
    b = (f(i-2)-f(i-1))*(f(i-2)-f(i));
    c = f(i-2)*f(i);
    d = (f(i-1)-f(i-2))*(f(i-1)-f(i));
    e = f(i-2)*f(i-1);
    f1 = (f(i) - f(i-2))*(f(i) - f(i-1));
    x(i+1) = a/b*x(i-2) + c/d*x(i-1) + e/f1*x(i);
    f(i+1) = F(x(i+1));
    Ea(i+1) = abs((x(i+1)-x(i))/x(i+1))*100;
    i=i+1;
end
x(i)