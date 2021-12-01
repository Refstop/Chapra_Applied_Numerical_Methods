clc;clear all;close all;
E = 200*10^9; I = 0.0003; w0 = 2.5*10^5; L = 3;
h = 0.125;
x = 0:h:L;
th = @(x) w0/(120*E*I*L)*(-5*x.^4+6*L^2*x.^2-L^4);
thf = th(x);
n = length(x);
%% (a) 보의 처짐
I = (x(n) - x(1)) * (thf(1)+2*sum(thf(2:n-1))+thf(n))/(2*n);
fprintf("(a) 보의 처짐: %.9f\n", I);

%% (b) 모멘트와 전단
M(1) = E*I*(thf(2)-thf(1))/h;
V(1) = E*I*(thf(3)-2*thf(2)+thf(1))/(h^2);
for xi = 2:n-1
    M(xi) = E*I*(thf(xi+1)-thf(xi-1))/(2*h);
    V(xi) = E*I*(thf(xi+1)-2*thf(xi)+thf(xi-1))/(h^2);
end
M(n) = E*I*(thf(n)-thf(n-1))/h;
V(n) = E*I*(thf(n)-2*thf(n-1)+thf(n-2))/(h^2);

fprintf("(b)\n");
for i = 1:length(x)
    fprintf("모멘트: %.9f\n", M(i));
    fprintf("전단력: %.9f\n\n", V(i));
end