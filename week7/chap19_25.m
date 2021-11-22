clc;clear all;close all;

z = [0 50 100 150 225 300 375 450 600];
F = [0 30 40 40 50 50 60 80 100];

%% (a) 빌딩에 가해지는 힘을 뉴턴 단위로 계산
xa1 = z(1:4); Fa1 = F(1:4);
n = length(xa1);
Ia1 = (xa1(n)-xa1(1))*(Fa1(1)+3*Fa1(2)+3*Fa1(3)+Fa1(4))/8;

xa2 = z(4:8); Fa2 = F(4:8);
n=length(xa2); h = (xa2(n)-xa2(1))/n;
odd = 0; even = 0;
for i = 2:n-1
    if rem(i,2) == 0
        even = even + Fa2(i);
    else
        odd = odd + Fa2(i);
    end
end
Ia2 = h/3*(Fa2(1)+4*odd+2*even+Fa2(n));

syms x;
xa3 = z(8:9); Fa3 = F(8:9);
Fa(x) = (Fa3(2) - Fa3(1)/(xa3(2)-xa3(1)))*(x-xa3(1))+Fa3(1);
Ia3 = int(Fa, xa3);
Ia = Ia1+Ia2+Ia3;
fprintf("힘의 크기(사다리꼴, Simpson 1/3 및 3/8 공식): %d\n", Ia);

%% (b) 힘의 작용선
fprintf("힘의 작용선:");
disp(z);