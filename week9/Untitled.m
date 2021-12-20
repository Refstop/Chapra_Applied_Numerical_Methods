clc;clear all;close all;

%% ODE 23-2
dydt = @(t,y) 10*exp(-(t-2)^2/(2*0.075^2))-0.6*y;

%% Euler
dt = 0.001; ti = 0; tf = 4;
t = [ti:dt:tf];
L = length(t);
if t(L)<tf
    t(L+1) = tf;
    L = L+1;
end
y0 = 0.5;
y(1) = y0;
for i = 1:L-1
    y(i+1) = y(i) + dydt(t(i), y(i))*dt;
end
figure; plot(t,y,'linewidth',5); grid on;
hold on;
dt = 0.1; ti = 0; tf = 4;
t = [ti:dt:tf];
L = length(t);
if t(L)<tf
    t(L+1) = tf;
    L = L+1;
end
y = 0;
y0 = 0.5;
y(1) = y0;
for i = 1:L-1
    y(i+1) = y(i) + dydt(t(i), y(i))*dt;
end
figure; plot(t,y,'ro'); grid on;

%% Heun
dt = 0.1; ti = 0; tf = 4;
t = [ti:dt:tf];
L = length(t);
if t(L)<tf
    t(L+1) = tf;
    L = L+1;
end
y = 0;
y0 = 0.5;
y(1) = y0;
Es = 0.1; 
for i = 1:L-1
    yo(i+1) = y(i) + dydt(t(i), y(i))*dt;
    Ea = 10; j=1;
    yc=yo(i+1);
    while Ea >= Es
        yc(j+1) = y(i) + ((dydt(t(i), y(i)) + dydt(t(i+1),yc(j)))/2)*dt;
        Ea = abs(yc(j+1)-yc(j));
        j=j+1;
    end
    y(i+1) = yc(j);
end
figure; plot(t,y); grid on;

%% RK 4
dt = 0.1; ti = 0; tf = 4;
t = [ti:dt:tf];
L = length(t);
if t(L)<tf
    t(L+1) = tf;
    L = L+1;
end
y = 0;
y0 = 0.5;
y(1) = y0;
for i = 1:L-1
    k1 = dydt(t(i), y(i));
    k2 = dydt(t(i)+(1/2)*dt, y(i)+(1/2)*k1*dt);
    k3 = dydt(t(i)+(1/2)*dt, y(i)+(1/2)*k2*dt);
    k4 = dydt(t(i)+dt, y(i)+k3*dt);
    y(i+1) = y(i)+1/6*(k1+2*k2+2*k3+k4)*dt;
end
figure; plot(t,y); grid on;

%% ODE23
dt = 0.1; ti = 0; tf = 4;

y = 0;
y0 = 0.5;

y(1) = y0;
DT(1) = dt;
Es = 0.001;
Tor = 0.00000000001;

E = 0;
T(1) = 0;
i = 1; 
while T <= 4
    j = 1;
    Next = 0;
    E = 0;
    while Next == 0
        k1 = dydt(T(i), y(i));
        k2 = dydt(T(i)+(1/2)*dt, y(i)+(1/2)*k1*dt);
        k3 = dydt(T(i)+(3/4)*dt, y(i)+(3/4)*k2*dt);
        k4 = dydt(T(i)+dt, y(i)+dt);
        y23 = y(i)+(1/9)*(2*k1+3*k2+4*k3)*dt;
        E(j+1) = abs(1/72*(-5*k1+6*k2+8*k3-9*k4)*dt);
        if j == 1
            Next = 0;
        else
            if E(j+1) > (Es+Tor) && E(j) > (Es+Tor)
                dt = dt*abs(Es/E(j+1))^0.2;
                Next = 0;
            elseif E(j+1) < (Es-Tor) && E(j) < (Es-Tor)
                dt = dt*abs(Es/E(j+1))^0.25;
                Next = 0;
            else
                Next = 1;
            end
        end
        j = j + 1;
    end
    y(i+1) = y23;
    T(i+1) = T(i) + dt;
    DT(i+1) = dt;
    i = i + 1;
end

figure(1); plot(T, y, 'ro'); grid on;
figure(2); plot(T, DT); grid on;