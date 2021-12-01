clc;clear all;close all;

f = @(x) exp(-2*x)-x;
xi = 2;

%% (a) 미적분학을 이용한 x=2에서의 도함수값
fd1 = @(x) -2*exp(-2*x)-1;
fprintf("(a) 미적분학을 이용한 x=2에서의 도함수값: %.6f\n", fd1(xi));

%% (b), (c), (d) 증분이 변하는 중심, 전향, 후향 유한차분 근사, 실제 해와의 비교
fprintf("(b), (c), (d)\n");
for dx = 0.5:-0.01:0.01
    xbc = [xi-dx xi xi+dx];
    fbc = f(xbc);
    cfd1 = (fbc(3)-fbc(1))/(2*dx);
    ffd1 = (-fbc(3)+4*fbc(2)-3*fbc(1))/(2*dx);
    bfd1 = (3*fbc(3)-4*fbc(2)+fbc(1))/(2*dx);
    fprintf("중심차분: %.6f\n", cfd1);
    fprintf("전향차분: %.6f\n", ffd1);
    fprintf("후향차분: %.6f\n", bfd1);
    fprintf("실제 해: %.6f\n\n", fd1(xi));
end

%% (d)