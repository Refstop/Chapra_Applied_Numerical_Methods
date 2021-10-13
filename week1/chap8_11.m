% 8-11
clc; clear all; close all;
k1 = 10; k2 = 40; k3 = 40; k4 = 10;
x1 = 0.05; x2 = 0.04; x3 = 0.03;
m1 = 1; m2 = 1; m3 = 1;
km_mat = [(k1+k2)/m1 -k2/m1 0;
    -k2/m2 (k2+k3)/m2 -k3/m2;
    0 -k3/m3 (k3+k4)/m3];
x = [x1; x2; x3];
X = -km_mat*x