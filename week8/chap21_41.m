clc;clear all;close all;
syms t;
v(t) = 2*t/sqrt(1+t^2);
a = diff(v,t);
tt = 5;
a_exact = a(tt);

h1 = 0.5; h2 = 0.25;
t1 = [tt-h1 tt tt+h1];
t2 = [tt-h2 tt tt+h2];
D05 = (v(t1(3))-v(t1(1)))/(2*h1);
D025 = (v(t2(3))-v(t2(1)))/(2*h2);
D = 4/3*D025-1/3*D05;
fprintf("물체의 가속도(엄밀): %.6f\n", a_exact);
fprintf("물체의 가속도(근사): %.6f\n", D);
e = abs((a_exact - D)/a_exact)*100;
fprintf("참 백분율 상대오차: %.6f%%\n", e);