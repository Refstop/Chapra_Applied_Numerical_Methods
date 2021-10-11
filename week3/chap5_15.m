clc; clear all;close all;
L = 6; E = 5*10^11; I = 3*10^(-4); w0 = 2.5*10^5;
F_coeffi = (w0/(120*E*I*L))*[-1 0 2*L^2 0 0 -L^4 0];
F_diff_coeffi = (w0/(120*E*I*L))*[-5 0 6*L^2 0 -L^4 0];
F = @(x) (w0/(120*E*I*L))*(-x.^5 + 2*L^2*x.^3 - L^4*x);
F_diff = @(x) (w0/(120*E*I*L))*(-5*x.^4 + 6*L^2*x.^2 - L^4);
xL = 0;
xU = 5.99;

%% incsearch
N = 10000;
x = linspace(xL, xU, N);
f = F_diff(x);
nb = 0; xb = [];
for k = 1:length(x)-1
    if sign(f(k)) ~= sign(f(k+1))
        nb = nb + 1;
        xb(nb) = (x(k)+x(k+1))/2;
    end
end
curve = F(xb);
fprintf("incsearch: %.16f\n", curve(find(max(abs(curve)))));

%% bisearch
Es = 0.001;
N = ceil(log2((xU - xL)/Es));

for i = 1:N
    xr(i) = (xU + xL)/2;
    SL = F_diff(xL)*F_diff(xr(i));
    if SL < 0
        xU = xr(i);
    else
        xL = xr(i);
    end
end

fprintf("bisearch: %.16f\n", max(F(xr(N))));

%% interpolation
Es = 0.001;
E = 1;
i = 0;
while(i<50)
    i = i + 1;
    xr(i) = xU - F_diff(xU)*(xL-xU)/(F_diff(xL)-F_diff(xU));
    SL = F_diff(xL)*F_diff(xr(i));
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

fprintf("interpolation: %.16f\n", max(F(xr(i))));
