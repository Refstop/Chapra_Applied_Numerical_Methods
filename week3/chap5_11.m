 clc;clear all;close all;
xLxU = [0 2;
    2 4.5;
    4.5 10;
    10 12];

for i = 1:4
    if i == 1
        M = @(x) -185*x;
    elseif i == 2
        M = @(x) 135*x+100;
    elseif i == 3
        M = @(x) -265*x+1450;
    else
        M = @(x) 100*x-1200;
    end
    xL = xLxU(i, 1);
    xU = xLxU(i, 2);
    %% incsearch
    N = 10000;
    x = linspace(xL, xU, N);
    f = M(x);
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
        SL = M(xL)*M(xr(i));
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
        xr(i) = xU - M(xU)*(xL-xU)/(M(xL)-M(xU));
        SL = M(xL)*M(xr(i));
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
    fprintf("\n");
end