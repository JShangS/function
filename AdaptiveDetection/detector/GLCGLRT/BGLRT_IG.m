clc
clear 
close all
%《Knowledge-based adaptive detection of radar targets 
%  in generalized Pareto clutter》
Na = 4;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
SNRout=0:1:20; % 输出SNR
cos2=0.9;
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
rou = 0.95;  %%协方差矩阵生成的迟滞因子
rouR = zeros(N,N);  %%真实的杂波协方差
L=round(2*N); 
theta_sig = 0.0;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
rou = 0.9:0.99;
for k = 1:length(rou)
    for i=1:N
        for j=1:N
            rouR(i,j,k)=rou(k)^abs(i-j)*exp(1j*2*pi*abs(i-j)*theta_sig);
        end
    end
end

