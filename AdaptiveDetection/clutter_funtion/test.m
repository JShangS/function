clc
clear 
close all
rou = 1;  %%协方差矩阵生成的迟滞因子
N = 8;
theta_sig = 0.1;
for i=1:N
    for j=1:N
        rouR(i,j)=rou^(i-j)*exp(1j*2*pi*(i-j)*theta_sig);%
    end
end
M = 2000;
ft = linspace(-0.5,0.5,M);
PSD=fun_PSD(rouR,ft);
plot(ft,abs(PSD))
