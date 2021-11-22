clc;clear all;close all;
t = [0 4 8 12 16 20 24 28 30];
v = [0 18 31 42 50 56 61 65 70];

%% (a) 30초까지 움직인 거리
ta1 = t(1:5); va1 = v(1:5);
n=length(ta1); h = (ta1(n)-ta1(1))/n;
odd = 0; even = 0;
for i = 2:n-1
    if rem(i,2) == 0
        even = even + va1(i);
    else
        odd = odd + va1(i);
    end
end
Ia1 = h/3*(va1(1)+4*odd+2*even+va1(n));

ta2 = t(5:8); va2 = v(5:8);
n = length(ta2);
Ia2 = (ta2(n)-ta2(1))*(va2(1)+3*va2(2)+3*va2(3)+va2(4))/8;

syms x;
ta3 = t(8:9); va3 = v(8:9);
va(x) = (va3(2) - va3(1)/(ta3(2)-ta3(1)))*(x-ta3(1))+va3(1);
Ia3 = int(va, ta3);
Ia = Ia1+Ia2+Ia3;
fprintf("움직인 거리(사다리꼴, Simpson 1/3 및 3/8 공식): %.6f\n", Ia);

%% (b) 평균 속도
avg_v = Ia/t(length(t));
fprintf("평균 속도: %.6f\n", avg_v);