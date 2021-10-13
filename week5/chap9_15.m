clc;clear all;close all;
Fbc = @(BD) -0.6*BD-74;
Fce = @(BC) 5/3*BC;
Fcd = @(CE) 0.8*CE+24;
Fbd = @(CD) -1.25*CD;
Fde = @(CE) -0.6*CE;
Fad = @(DE, BD) DE - 0.6*BD;
Fey = @(CE) -0.8*CE;
Fab = @(BD) -0.8*BD;
Fay = @(AB) -AB;
Fax = @(AD) -AD;

A = [0 0 1 0 0 0 0 1 0 0;
     1 0 0 0 0 0 0 0 1 0;
     0 1 0 0.6 0 0 0 0 0 0;
     -1 0 0 -0.8 0 0 0 0 0 0;
     0 -1 0 0 0 0 0.6 0 0 0;
     0 0 0 0 -1 0 -0.8 0 0 0;
     0 0 -1 -0.6 0 1 0 0 0 0;
     0 0 0 0.8 1 0 0 0 0 0;
     0 0 0 0 0 -1 -0.6 0 0 0;
     0 0 0 0 0 0 0.8 0 0 1];
b = [0; 0; -74; 0; 0; 24; 0; 0; 0; 0];
X = inv(A)*b

%% 연속대입법(Gauss Seidal)
% AB = 0; BC = 0; AD = 0; BD = 0; CD = 0; DE = 0; CE = 0; Ax = 0; Ay = 0; Ey = 0;
BD = 0;
es = 0.001;
while(1)
    BDold = BD;
    BC = Fbc(BD);
    CE = Fce(BC);
    CD = Fcd(CE);
    BD = Fbd(CD);
    DE = Fde(CE);
    AD = Fad(DE, BD);
    Ey = Fey(CE);
    AB = Fab(BD);
    Ay = Fay(AB);
    Ax = Fax(AD);
    ea = abs((BD-BDold)/BD)*100;
    if ea <= es
        break;
    end
end
fprintf("연속대입법(Gauss Seidal):\n");
fprintf("AB: %.8f  BC: %.8f  AD: %.8f BD: %.8f CD: %.8f\n", AB, BC, AD, BD, CD);
fprintf("DE: %.8f  CE: %.8f  Ax: %.8f Ay: %.8f Ey: %.8f\n", DE, CE, Ax, Ay, Ey);

%% 연속대입법(Jacobi)
AB = 0; BC = 0; AD = 0; BD = 0; CD = 0; DE = 0; CE = 0; Ax = 0; Ay = 0; Ey = 0;
es = 0.001;
while(1)
    ABold = AB; BCold = BC; ADold = AD; BDold = BD; CDold = CD;
    DEold = DE; CEold = CE; Axold = Ax; Ayold = Ay; Eyold = Ey;
    BC = Fbc(BDold);
    CE = Fce(BCold);
    CD = Fcd(CEold);
    BD = Fbd(CDold);
    DE = Fde(CEold);
    AD = Fad(DEold, BDold);
    Ey = Fey(CEold);
    AB = Fab(BDold);
    Ay = Fay(ABold);
    Ax = Fax(ADold);
    ea = abs((BD-BDold)/BD)*100;
    if ea <= es
        break;
    end
end
fprintf("연속대입법(Jacobi):\n");
fprintf("AB: %.8f  BC: %.8f  AD: %.8f BD: %.8f CD: %.8f\n", AB, BC, AD, BD, CD);
fprintf("DE: %.8f  CE: %.8f  Ax: %.8f Ay: %.8f Ey: %.8f\n", DE, CE, Ax, Ay, Ey);