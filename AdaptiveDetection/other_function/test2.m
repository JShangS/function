clc
clear 
close all
% rou = 1;  %%协方差矩阵生成的迟滞因子
% N = 8;
% theta_sig = 0.1;
% for i=1:N
%     for j=1:N
%         rouR(i,j)=rou^(i-j)*exp(1j*2*pi*(i-j)*theta_sig);%
%     end
% end
n = 4; %几倍的样本
str_train = 'g';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 3; %%%IG的选项，1为每个距离单元IG纹理都不同
sigma_t = 1;
rou = 0.95;  %%协方差矩阵生成的迟滞因子
%%%Pd_CLGLRT_2Kmu1lambda3s0.1o1_p：2K：训练单元数目，mu，lambda，s：失配向量方差，
%%o1:opt=1，p：IG纹理复合高斯
%%%%假设参数设置
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
L=round(n*N); 
for i=1:N
    for j=1:N
        rouR(i,j)=rou^abs(i-j);%*exp(1j*2*pi*abs(i-j)*theta_sig);
    end
end
Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
R_SCM = (fun_SCMN(Train));
R_NSCM = fun_NSCMN(Train);
R_AML = fun_AML(Train);
M = 1000;
ft = linspace(-0.5,0.5,M);
PSD=fun_PSD(rouR,ft);
% PSD=PSD/max(abs(PSD));
% PSD=fun_value2dB(PSD);
plot(ft,(PSD))
hold on
plot(ft,abs(fun_PSD(R_SCM,ft)),'r')
plot(ft,(fun_PSD(R_NSCM,ft)),'k.-')
plot(ft,(fun_PSD(R_AML,ft)),'g.-')
