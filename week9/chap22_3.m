clc;clear all;close all;
syms y(t);
dydt = diff(y,t) == -y+t^2;
ti = 0; tf = 3; dt = 0.001;
t = [ti:dt:tf]; y0 = 1;

%% 수치 미방해석 공통
dydt = @(t,y) -y+t^2;
yy=0; yy(1) = y0;

%% (a) 반복이 없는 Heun법
for i = 1:length(t)-1
    k1 = dydt(t(i), yy(i));
    k2 = dydt(t(i)+dt, yy(i)+k1*dt);
    yy(i+1) = yy(i) + (0.5*k1+0.5*k2)*dt;
end
figure;
plot(t, yy, 'go'); grid on;
hold on;
%% (b) 조건부 Heun법
Es = 0.001;
for i = 1:length(t)-1
    yo(i+1) = yy(i) + dydt(t(i), yy(i))*dt;
    Ea = 10; j = 1;
    yc=yo(i+1);
    while Ea >= Es
        yc(j+1) = yy(i) + (dydt(t(i),yy(i)) + dydt(t(i+1),yc(j)))/2*dt;
        Ea = abs((yc(j+1)-yc(j))/yc(j+1))*100;
        j = j + 1;
    end
    yy(i+1) = yc(j);
end
plot(t, yy); grid on;
%% (c) 중점법
for i = 1:length(t)-1
    k1 = dydt(t(i), yy(i));
    k2 = dydt(t(i)+(1/2)*dt, yy(i)+(1/2)*k1*dt);
    yy(i+1) = yy(i) + k2*dt;
end
% plot(t, yy, 'ko'); grid on;
%% (d) Ralston법
for i = 1:length(t)-1
    k1 = dydt(t(i), yy(i));
    k2 = dydt(t(i)+(2/3)*dt, yy(i)+(2/3)*k1*dt);
    yy(i+1) = yy(i) + (0.25*k1+0.75*k2)*dt;
end
% plot(t, yy); grid on;