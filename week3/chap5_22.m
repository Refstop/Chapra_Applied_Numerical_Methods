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
Es = 0.001;
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
        E = abs((xr(i)-xr(i-1))/xr(i))*100;
        if E < Es
            break; 
        end
    end
    
end
fprintf("interpolation: %.16f\n", xr(i));
