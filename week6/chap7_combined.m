clc;clear all;close all;
% F = @(x) x^2/10-2*sin(x);
F = @(x) 4*x - 1.8*x^2 + 1.2*x^3 - 0.3*x^4;

%% Brent(2차 보간법 + 황금분할탐색)
i = 0;
x1 = 0;
x2 = 1;
x3 = 4;

Es = 0.001; Ea = 1;

while(Ea > Es)
    term1 = (x2-x1)*(F(x2)-F(x3));
    term2 = (x2-x1)*term1;
    term3 = (x2-x3)*(F(x2)-F(x1));
    term4 = (x2-x3)*term3;
    xopt_2 = x2-0.5*(term2-term4)/(term1-term3);
    if xopt_2 < x2
        Ea = abs((x2-x1)/xopt_2)*100;
        x3 = x2; x2 = xopt_2;
    else
        Ea = abs((x3-x2)/xopt_2)*100;
        x1 = x2; x2 = xopt_2;
    end
    
    if i>=100
        break;
    end
    i=i+1;
end
fprintf('(1) 2차 보간법 종료: %.8f\n', xopt_2);

dx = 0.01;
xl = xopt_2-dx; xu = xopt_2+dx;
Es = 0.001; Ea = 1;
phi = (1+sqrt(5))/2;
i = 0;
d = xu-xl;
while(Ea > Es)
    d = (phi-1)*d;
    x1 = xl+d;
    x2 = xu-d;
    if F(x1) < F(x2)
        xopt_g = x1; xl = x2;
        Ea = max(abs((x1-x2)/xopt_g)*100, abs((xu-x1)/xopt_g)*100);
    else
        xopt_g = x2; xu = x1;
        Ea = max(abs((x1-x2)/xopt_g)*100, abs((x2-xl)/xopt_g)*100);
    end
    if i>=100
        break;
    end
    i=i+1;
end
fprintf('(2) 황금분할탐색: %.8f\n', xopt_g);