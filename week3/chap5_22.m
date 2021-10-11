clc; clear all;close all;
g = 9.81; u = 1800; m0 = 160000; q = 2600; v = 750;
F = @(t) u*log(m0 ./ (m0 - q*t))-g*t-v;
xL = 10;
xU = 50;

%% incsearch
N = 10000;
x = linspace(xL, xU, N);
f = F(x);
nb = 0; xb = [];
for k = 1:length(x)-1
    if sign(f(k)) ~= sign(f(k+1))
        nb = nb + 1;
        xb(nb) = (x(k)+x(k+1))/2;
    end
end
fprintf("incsearch: %.16f\n", xb);

%% bisearch
Es = 0.01;
N = ceil(log2((xU - xL)/Es));

for i = 1:N
    xr(i) = (xU + xL)/2;
    SL = F(xL)*F(xr(i));
    if SL < 0
        xU = xr(i);
    else
        xL = xr(i);
    end
end
fprintf("bisearch: %.16f\n", xr(N));

%% interpolation
Es = 0.01;
E = 1;
i = 0;
while(i<50)
    i = i + 1;
    xr(i) = xU - F(xU)*(xL-xU)/(F(xL)-F(xU));
    SL = F(xL)*F(xr(i));
    if SL < 0
        xU = xr(i);
    else
        xL = xr(i);
    end
    if i > 1
        E = abs((xr(i)-xr(i-1))/xr(i));
        if E < Es
            break; 
        end
    end
    
end
fprintf("interpolation: %.16f\n", xr(i));
