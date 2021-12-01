clc;clear all;close all;
t = [0 2 4 6 8 10 12 14 16];
x = [0 0.7 1.8 3.4 5.1 6.3 7.3 8.0 8.4];

xi = 6;
h = 2;
%% (a) 중심차분
va = (x(xi+1)-x(xi-1))/(2*h);
aa = (x(xi+1)-2*x(xi)+x(xi-1))/(h^2);
fprintf("(a)\n");
fprintf("중심차분 속도 v=%.6f\n", va);
fprintf("중심차분 가속도 a=%.6f\n", aa);

%% (b) 전향차분
vb = (-x(xi+2)+4*x(xi+1)-3*x(xi))/(2*h);
ab = (-x(xi+3)+4*x(xi+2)-5*x(xi+1)+2*x(xi))/(h^2);
fprintf("(b)\n");
fprintf("전향차분 속도 v=%.6f\n", vb);
fprintf("전향차분 가속도 a=%.6f\n", ab);

%% (c) 후향차분
vc = (3*x(xi)-4*x(xi-1)+x(xi-2))/(2*h);
ac = (2*x(xi)-5*x(xi-1)+4*x(xi-2)-x(xi-3))/(h^2);
fprintf("(c)\n");
fprintf("후향차분 속도 v=%.6f\n", vc);
fprintf("후향차분 가속도 a=%.6f\n", ac);