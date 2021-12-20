clc;clear all;close all;
syms x
F(x) = 8+4*cos(x);
a = 0; b = pi/2;

%% (a) 해석적인 방법
Ia = int(F, [a b]);
fprintf("해석적인 방법: %.6f\n", Ia);

%% (b) 단일 구간에 대한 사다리꼴 공식
Ib = (b-a)*(F(b)+F(a))/2;
e = abs((Ia - Ib)/Ia)*100;
fprintf("단일 구간에 대한 사다리꼴 공식: %.6f    오차율: %.3f%%\n", Ib, e);

%% (c) 합성 사다리꼴 공식
n = 2;
xc = linspace(a, b, n);
Ic = (xc(n) - xc(1)) * (F(xc(1))+2*sum(F(xc(2:n-1)))+F(xc(n)))/(2*n);
e = abs((Ia - Ic)/Ia)*100;
fprintf("합성 사다리꼴 공식(n=2): %.6f    오차율: %.3f%%\n", Ic, e);
n = 4;
xc = linspace(a, b, n);
Ic = (xc(n) - xc(1)) * (F(xc(1))+2*sum(F(xc(2:n-1)))+F(xc(n)))/(2*n);
e = abs((Ia - Ic)/Ia)*100;
fprintf("합성 사다리꼴 공식(n=4): %.6f    오차율: %.3f%%\n", Ic, e);

%% (d) 단일 구간에 대한 Simpson 1/3 공식
x1 = (a+b)/2;
Id = (b-a)*(F(a)+4*F(x1)+F(b))/6;
e = abs((Ia - Id)/Ia)*100;
fprintf("단일 구간에 대한 Simpson 1/3 공식: %.6f    오차율: %.3f%%\n", Id, e);

%% (e) 합성 Simpson 1/3 공식
n=4; xe = linspace(a, b, n); h = (b-a)/n;
odd = 0; even = 0;
for i = 2:n-1
    if rem(i,2) == 0
        even = even + F(xe(i));
    else
        odd = odd + F(xe(i));
    end
end
Ie = h/3*(F(xe(1))+4*odd+2*even+F(xe(n)));
e = abs((Ia - Ie)/Ia)*100;
fprintf("합성 Simpson 1/3 공식: %.6f    오차율: %.3f%%\n", Ie, e);

%% (f) Simpson 3/8 공식
xf = linspace(a, b, 4);
If = (b-a)*(F(xf(1))+3*F(xf(2))+3*F(xf(3))+F(xf(4)))/8;
e = abs((Ia - If)/Ia)*100;
fprintf("Simpson 3/8 공식: %.6f    오차율: %.3f%%\n", If, e);

%% (g) 합성 Simpson 공식
n=5; xg = linspace(a, b, n); h = (b-a)/n;
odd = 0; even = 0;
for i = 2:n-1
    if rem(i,2) == 0
        even = even + F(xg(i));
    else
        odd = odd + F(xg(i));
    end
end
Ig = h/3*(F(xg(1))+4*odd+2*even+F(xg(n)));
e = abs((Ia - Ig)/Ia)*100;
fprintf("합성 Simpson 공식: %.6f    오차율: %.3f%%\n", Ig, e);