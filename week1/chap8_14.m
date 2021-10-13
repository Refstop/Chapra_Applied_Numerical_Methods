% 8-14
clc; clear all; close all;
R12 = 10; R23 = 20; R34 = 2; R45 = 5; R56 = 25; R25 = 5;
V1 = 150; V6 = 0;
% i12, i23(=i34=i45), i25, i56
% clockwise, counter-clockwise
A = [0 R23+R34+R45 R25 0;
    R12 0 R25 R56;
    1 1 -1 0;
    0 1 -1 1];
b = [0; -(V6-V1); 0; 0];
I = inv(A)*b;
i12 = I(1)
i23 = I(2)
i34 = I(2)
i45 = I(2)
i25 = I(3)
i56 = I(4)