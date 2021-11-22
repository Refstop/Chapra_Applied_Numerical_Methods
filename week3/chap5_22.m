clc; clear all;close all;
g = 9.81; u = 1800; m0 = 160000; q = 2600; v = 750;
F = @(t) u*log(m0 ./ (m0 - q*t))-g*t-v;
xL = 10;
xU = 50;

%% incsearch
N = 10000;
x = linspace(xL, xU, N);
f = F(x);
xb = [];
for k = 1:length(x)-1
    if sign(f(k)) ~= sign(f(k+1))
        xb = (x(k)+x(k+1))/2;
    end
end
fprintf("incsearch: %.16f\n", xb);

%% bisearch
Es = 0.001;
N = ceil(log2((xU - xL)/Es));

for i = 1:N
    xr(i) = (xU + xL)/2;
    SL = F(xL)*F(xr(i));
    if SL < 0
        xU = xr(i);
    elseif SL > 0
        xL = xr(i);
    else
        break;
    end
end
fprintf("bisearch: %.16f\n", xr(N));

%% interpolation
xL = 10; xU = 50;
xr = 0;
i = 0; es = 0.001; ea = 1;

while(ea>es)
    xrold = xr;
    xr = xU - F(xU)*(xL-xU)/(F(xL)-F(xU));
    SL = F(xL)*F(xr);
    if SL < 0
        xU = xr;
    elseif SL > 0
        xL = xr;
    else
        break;
    end
    ea = abs((xr-xrold)/xr)*100;
    if i >= 100
        break;
    end
    i=i+1;
end
fprintf("interpolation: %.16f\n", xr);
