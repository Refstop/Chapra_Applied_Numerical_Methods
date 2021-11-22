clc;clear all;close all;
Fbc = @(BD) -0.6*BD-74;
Fce = @(BC) 5/3*BC;
Fcd = @(CE) -0.8*CE-24;
Fbd = @(CD) -5/4*CD;
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
BD = 1;
es = 0.001; i=0;
while(i<=1000)
    BDold = BD;
    BC = Fbc(BD);
    CE = Fce(BC);
    CD = Fcd(CE);
    BD = Fbd(CD);
    DE = Fde(CE);
    AD = Fad(DE, BD);
    Ax = Fax(AD);
    Ey = Fey(CE);
    AB = Fab(BD);
    Ay = Fay(AB);
    
    ea = abs((BD-BDold)/BD)*100;
    if ea <= es
        break;
    end
    i=i+1;
end
if i>=1000
    fprintf("방정식의 기울기 값의 절댓값이 1보다 큰 경우가 존재하므로 발산한다.\n");
else
    fprintf("연속대입법(Gauss Seidal):\n");
    fprintf("AB: %.8f  BC: %.8f  AD: %.8f BD: %.8f CD: %.8f\n", AB, BC, AD, BD, CD);
    fprintf("DE: %.8f  CE: %.8f  Ax: %.8f Ay: %.8f Ey: %.8f\n", DE, CE, Ax, Ay, Ey);
end

%% 연속대입법(Jacobi)
AB = 1; BC = 1; AD = 1; BD = 1; CD = 1; DE = 1; CE = 1; Ax = 1; Ay = 1; Ey = 1;
es = 0.001; i=0;
while(i<=1000)
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
    i=i+1;
end
if i>=1000
    fprintf("방정식의 기울기 값의 절댓값이 1보다 큰 경우가 존재하므로 발산한다.\n");
else
    fprintf("연속대입법(Jacobi):\n");
    fprintf("AB: %.8f  BC: %.8f  AD: %.8f BD: %.8f CD: %.8f\n", AB, BC, AD, BD, CD);
    fprintf("DE: %.8f  CE: %.8f  Ax: %.8f Ay: %.8f Ey: %.8f\n", DE, CE, Ax, Ay, Ey);
end
