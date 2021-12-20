clc;clear all;close all;

%% euler 2계 미분방정식
ti = 0; tf = 4; dt = 0.1;
yinit = [1 0];

t = [ti:dt:tf];
n = length(t);
y = [yinit;zeros(n-1,2)];

for i = 1:n-1
    y(i+1,:) = y(i,:) + func(t(i),y(i,:))*dt;
end
figure;
plot(t,y(:,1)); grid on;
hold on;

%% RK 4차 2계 미분방정식
ti = 0; tf = 4; dt = 0.1;
yinit = [1 0];

t = [ti:dt:tf];
n = length(t);
y = [yinit;zeros(n-1,2)];

k=zeros(4,2);

for i = 1:n-1
    k(1,:)=func(t(i),y(i,:));
    k(2,:)=func(t(i)+dt/2,y(i,:)+(dt/2)*k(1,:));
    k(3,:)=func(t(i)+dt/2,y(i,:)+(dt/2)*k(2,:));
    k(4,:)=func(t(i)+dt,y(i,:)+dt*k(3,:));
    y(i+1,:) = y(i,:) + (dt/6)*(k(1,:)+2*k(2,:)+2*k(3,:)+k(4,:));
end
plot(t,y(:,1)); grid on;

%% y=cos 3t
ti = 0; tf = 4; dt = 0.1;
t = [ti:dt:tf];

y = cos(3*t);
plot(t,y,'o'); grid on;

function dy = func(t, y)
  dy = [y(2),-9*y(1)];
end
