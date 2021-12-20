clc;clear all;close all;
syms y(x);
dydx = diff(y,x) == (1+2*x)*sqrt(y);
ti = 0; tf = 1; dx = 0.01;
xx = [ti:dx:tf]; y0 = 1;

%% (a) 해석적인 방법
cond = y(0) == 1;
ySol(x) = dsolve(dydx, cond);
yy = ySol(xx);
yy1 = double(yy{1});
yy2 = double(yy{2});
figure; plot(xx, yy2); grid on;
hold on;

%% 수치 미방해석 공통
dydx = @(x,y) (1+2*x)*sqrt(y);
yy=0; yy(1) = y0;

%% (b) Euler 방법
for i = 1:length(xx)-1
    yy(i+1) = yy(i) + dydx(xx(i), yy(i))*dx;
end
plot(xx, yy); grid on;
%% (c) 반복이 없는 Heun법
for i = 1:length(xx)-1
    k1 = dydx(xx(i), yy(i));
    k2 = dydx(xx(i)+dx, yy(i)+k1*dx);
    yy(i+1) = yy(i) + (0.5*k1+0.5*k2)*dx;
end
plot(xx, yy); grid on;
%% (d) Ralston법
for i = 1:length(xx)-1
    k1 = dydx(xx(i), yy(i));
    k2 = dydx(xx(i)+(2/3)*dx, yy(i)+(2/3)*k1*dx);
    yy(i+1) = yy(i) + (0.25*k1+0.75*k2)*dx;
end
plot(xx, yy); grid on;
%% (e) 4차 RK법
for i = 1:length(xx)-1
    k1 = dydx(xx(i), yy(i));
    k2 = dydx(xx(i)+0.5*dx, yy(i)+0.5*k1*dx);
    k3 = dydx(xx(i)+0.5*dx, yy(i)+0.5*k2*dx);
    k4 = dydx(xx(i)+dx, yy(i)+k3*dx);
    yy(i+1) = yy(i) + (1/6)*(k1+2*k2+2*k3+k4)*dx;
end
plot(xx, yy); grid on;