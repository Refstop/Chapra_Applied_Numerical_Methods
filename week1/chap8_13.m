% 8-13
clc; clear all; close all;
g = 9.81; k = 10;
m1 = 2; m2 = 3; m3 = 2.5;

A = [k 0; k -k];
b = [-m1*g; m2*g];
dx = inv(A)*b