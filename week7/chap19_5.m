clc;clear all;close all;
syms x
F(x) = exp(-x);
a = 0; b = 1.2;
xx = [0 0.1 0.3 0.5 0.7 0.95 1.2];

%% (a) 해석적인 방법
Ia = int(F, [a b]);
fprintf("해석적인 방법: %.6f\n", Ia);

%% (b) 사다리꼴 공식
Fb(x) = (F(b) - F(a)/(b-a))*(x-a)+F(a);
Ib = int(Fb, [a b]);
e = abs((Ia - Ib)/Ia)*100;
fprintf("단일 구간에 대한 사다리꼴 공식: %.6f    오차율: %.3f%%\n", Ib, e);

%% (c) 사다리꼴 공식과 Simpson 공식의 조합
xc = xx(1:5);
n = length(xc);
Ic1 = (xc(n) - xc(1)) * (F(xc(1))+2*sum(F(xc(2:n-1)))+F(xc(n)))/(2*n);
a = xx(5); b = xx(7);
x1 = (a+b)/2;
Ic2 = (b-a)*(F(a)+4*F(x1)+F(b))/6;
Ic = Ic1 + Ic2;
e = abs((Ia - Ic)/Ia)*100;
fprintf("사다리꼴 공식과 Simpson 공식의 조합: %.6f    오차율: %.3f%%\n", Ic, e);