clc;clear all;close all;
x = [-2 -1.5 -1 -0.5 0 0.5 1 1.5 2];
f = [0.05399 0.12952 0.24197 0.35207 0.39894 0.35207 0.24197 0.12952 0.05399];

h = 0.5;
fd1old = 0; x_infle = 0;

for xi = 1:length(x)-1
    fd1 = (f(xi+1)-f(xi))/h;
    if fd1*fd1old < 0
        x_infle = xi;
        break;
    end
    fd1old = fd1;
end
fprintf("변곡점 = %.6f\n", x(x_infle));