clc;clear all;close all;
Fx = @(y) (y+1)^0.5;
Fy = @(x) (5-x^2)^0.5;

%% 연속대입법
xr = 1.5; yr = 1.5;
es = 0.001;
while(1)
    xrold = xr;
    yrold = yr;
    xr = Fx(yr);
    yr = Fy(xr);
    ea = abs((xr-xrold)/xr)*100;
    if ea <= es
        break;
    end
end
fprintf("연속대입법(Gauss Seidal): %.8f %.8f\n", xr, yr);

%% Newton-Raphson
F1 = @(x,y) x^2+y^2-5;
F2 = @(x,y) x^2-y-1;
dF1dx = @(x) 2*x;
dF1dy = @(y) 2*y;
dF2dx = @(x) 2*x;
dF2dy = @(x) -1;

xyr = [1.5; 1.5];
es=0.001;
while(1)
    xrold = xyr(1);
    J = [dF1dx(xyr(1)) dF1dy(xyr(2)); dF2dx(xyr(1)) dF2dy(xyr(2))];
    F = [F1(xyr(1),xyr(2)); F2(xyr(1),xyr(2))];
    dx = J\F;
    xyr = xyr - dx;
    ea = abs((xyr(1)-xrold)/xyr(1))*100;
    if ea<=es
        break;
    end
end
fprintf("Newton Raphson: %.8f %.8f\n", xyr(1), xyr(2));