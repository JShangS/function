%%%几种协方差估计的比较
clc
clear 
close all
Na = 4;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
rou = 0.95;  %%协方差矩阵生成系数
rouR = zeros(N,N);
L=round(2*N); 
% theta_sig = 0.1;
% nn = 0:N-1;
% s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
for i=1:N
    for j=1:N
        rouR(i,j)=rou^abs(i-j);
    end
end
t = 0.1*rand(N,1);
R_KA = rouR.*(1+t*t');
irouR=inv(rouR);
rouR_half=rouR^0.5;
MC = 1000;
ALL = sum(sum(abs(rouR)));
error_SCM = zeros(MC,1);
error_NSCM = zeros(MC,1);
error_AML = zeros(MC,1);
error_AIWCM = zeros(MC,1);
error_CC_SCM = zeros(MC,1);
error_CC_NSCM = zeros(MC,1);
error_CC_AML = zeros(MC,1);
for i = 1:MC
    i
    %%产生杂波和噪声
    X = (randn(N,L)+1i*randn(N,L))/sqrt(2);  % 产生方差为1的复高斯白噪声 % Rwhite1=1/snapshot1*X1*X1'; eig(Rwhite1); % round(mean(abs(eig(Rwhite1)))) == 1
    Train = rouR_half*X;%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    W=(randn(N,1)+1i*randn(N,1))/sqrt(2); % 1i == -i
    x0=rouR_half*W;%+pp; % noise=(randn(N,1)+j*randn(N,1))/sqrt(2);  % 接收信号仅包括杂波和噪声
    %%协方差估计
    %%SCM
    R_SCM = fun_SCM(Train);
    %%NSCM
    R_NSCM = fun_NSCM(Train);
    %%AML
    R_AML = fun_AML(Train);
    %%AIWCM
%     R_AIWCM = fun_AIWCM(Train,x0);
    %%CC+R_SCM
    R_CC_SCM = fun_CC(Train,R_SCM,R_KA);
    %%CC+R_NSCM
    R_CC_NSCM = fun_CC(Train,R_NSCM,R_KA);
    %%CC+R_AML
    R_CC_AML = fun_CC(Train,R_AML,R_KA);
    %%%误差比较
    error_SCM(i) = sum(sum(abs(R_SCM-rouR)))/ALL;
    error_NSCM(i) = sum(sum(abs(R_NSCM-rouR)))/ALL;
    error_AML(i) = sum(sum(abs(R_AML-rouR)))/ALL;
%     error_AIWCM = error_AIWCM + sum(sum(abs(R_AIWCM-rouR)))/ALL/MC;
    error_CC_SCM(i) = sum(sum(abs(R_CC_SCM-rouR)))/ALL;
    error_CC_NSCM(i) = sum(sum(abs(R_CC_NSCM-rouR)))/ALL;
    error_CC_AML(i) = sum(sum(abs(R_CC_AML-rouR)))/ALL;
end
mean_error_SCM = mean(error_SCM);
mean_error_NSCM = mean(error_NSCM);
mean_error_AML = mean(error_AML);
mean_error_CC_SCM = mean(error_CC_SCM);
mean_error_CC_NSCM = mean(error_CC_NSCM);
mean_error_CC_AML = mean(error_CC_AML);

var_error_SCM = var(error_SCM);
var_error_NSCM = var(error_NSCM);
var_error_AML = var(error_AML);
var_error_CC_SCM = var(error_CC_SCM);
var_error_CC_NSCM = var(error_CC_NSCM);
var_error_CC_AML = var(error_CC_AML);

figure(1)
subplot(3,2,1)
plot(error_CC_NSCM,'g')
legend('CCNSCM')
axis([0,MC,0,1])
subplot(3,2,2)
plot(error_SCM,'r')
legend('SCM')
axis([0,MC,0,1])
subplot(3,2,3)
plot(error_NSCM,'b')
legend('NSCM')
axis([0,MC,0,1])
subplot(3,2,4)
plot(error_AML,'k')
legend('AML')
axis([0,MC,0,1])
subplot(3,2,5)
plot(error_CC_SCM,'y')
legend('CCSCM')
axis([0,MC,0,1])
subplot(3,2,6)
plot(error_CC_AML,'c')
legend('CCAML')
axis([0,MC,0,1])




