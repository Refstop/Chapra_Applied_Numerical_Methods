clc;clear all;close all;
Fx1 = @(x2, x3) (27 - 2*x2 + x3) / 10;
Fx2 = @(x1, x3) (61.5 - 3*x1 + 2*x3) / 5;
Fx3 = @(x1, x2) -(21.5 + x1 + x2) / 6;

A = [10 2 -1;
     -3 -5 2;
      1 1 6];
b = [27; -61.5; -21.5];
X = inv(A)*b;
fprintf("역함수: %.8f %.8f %.8f\n", X(1), X(2), X(3));

%% 연속대입법(Gauss Seidal)
xr1 = 0; xr2 = 0; xr3 = 0;
es = 0.001;
while(1)
    xr1old = xr1;
    xr1 = Fx1(xr2, xr3);
    xr2 = Fx2(xr1, xr3);
    xr3 = Fx3(xr1, xr2);
    ea = abs((xr1-xr1old)/xr1)*100;
    if ea <= es
        break;
    end
end
fprintf("연속대입법(Gauss Seidal): %.8f %.8f %.8f\n", xr1, xr2, xr3);

%% 연속대입법(Jacobi)
xr1 = 0; xr2 = 0; xr3 = 0;
es = 0.001;
while(1)
    xr1old = xr1;
    xr2old = xr2;
    xr3old = xr3;
    xr1 = Fx1(xr2old, xr3old);
    xr2 = Fx2(xr1old, xr3old);
    xr3 = Fx3(xr1old, xr2old);
    ea = abs((xr1-xr1old)/xr1)*100;
    if ea <= es
        break;
    end
end
fprintf("연속대입법(Jacobi): %.8f %.8f %.8f\n", xr1, xr2, xr3);
