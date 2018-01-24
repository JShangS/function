%%%协方差估计的信息几何方法，采用MAM模型为协方差的先验协方差。
clc
clear 
close all
rou1=0.95;
rou2=0.940;
N=8;
L=round(1*N);
LL=1000;
R1=fun_rho(rou1,N);
R2=fun_rho(rou2,N);
distance1x0=0;
distance2x0=0;
for i =1:LL
    i
    X = fun_TrainData_gauss(N,L,R1);
    R_x0=fun_SCMN(X);
    distance1x0 = distance1x0+fun_ReimanDistance(R1,R_x0)/LL;
    distance2x0 = distance2x0+fun_ReimanDistance(R2,R_x0)/LL;
end
ss=exp(-distance1x0)+exp(-distance2x0);
ratio1 = exp(-distance1x0)/ss;
ratio2 = exp(-distance2x0)/ss;