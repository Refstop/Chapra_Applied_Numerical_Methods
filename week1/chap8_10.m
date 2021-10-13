% 8-10
clc; clear all; close all;
% F1 F2 F3 H2 V2 V3
A = [cosd(30) 0 -cosd(60) 0 0 0;
    sind(30) 0 sind(60) 0 0 0;
    cosd(30) 1 0 1 0 0;
    sind(30) 0 0 0 1 0;
    0 1 cosd(60) 0 0 0;
    0 0 sind(60) 0 0 1];
b = [1000; 1000; 0; 0; 0; 0];
X = inv(A)*b;
F1 = X(1)
F2 = X(2)
F3 = X(3)
H2 = X(4)
V2 = X(5)
V3 = X(6)