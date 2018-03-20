%%%几种协方差估计的比较
clc
clear 
close all
Na = 4;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
rou = 0.95;  %%协方差矩阵生成的迟滞因子
rouR = zeros(N,N);  %%真实的杂波协方差
L=round(1.5*N); 
% theta_sig = 0.1;
% nn = 0:N-1;
% s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
for i=1:N
    for j=1:N
        rouR(i,j)=rou^abs(i-j);
    end
end
irouR=inv(rouR);
t = normrnd(1,0.01,N,1);%%0~0.5%%失配向量
R_KA = rouR.*(t*t');
rouR_half=rouR^0.5;
MC = 1000;
ALL = sum(sum(sqrt(abs(rouR).^2)));%真正协方差的范数
error_SCM = zeros(MC,1);
error_NSCM = zeros(MC,1);
error_AML = zeros(MC,1);
error_AIWCM = zeros(MC,1);
error_CC_SCM = zeros(MC,1);
error_CC_NSCM = zeros(MC,1);
error_CC_AML = zeros(MC,1);
error_CC_AIWCM = zeros(MC,1);
error_KCC_AIWCM = zeros(MC,1);
h = waitbar(0,'Please wait...');
for i = 1:MC
    waitbar(i/MC,h,sprintf([num2str(i/MC*100),'%%']));
    %%产生杂波和噪声
    Train = fun_TrainData('p',N,L,rouR,3,1,1);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = fun_TrainData('p',N,1,rouR,3,1,1); % 接收信号仅包括杂波和噪声
    Sx0 = (x0*x0');
%%协方差估计
    %%SCM，采样协方差
    R_SCM = abs(fun_SCM(Train));
    %%NSCM，归一化采样协方差
    R_NSCM = abs(fun_NSCM(Train));
    %%AML，近似最大似然
    R_AML = abs(fun_AML(Train));
    %%AIWCM，自适应迭代协方差估计
    R_AIWCM =abs(fun_AIWCM(Train,x0));
    %%CC+R_SCM，凸优化的知识辅助+协方差估计算法的训练数据协方差
    R_CC_SCM = abs(fun_CC(Train,R_SCM,R_KA));
    %%CC+R_NSCM
    R_CC_NSCM = abs(fun_CC(Train,R_NSCM,R_KA));
    %%CC+R_AML
    R_CC_AML = abs(fun_CC(Train,R_AML,R_KA));
    %%CC+R_AML
    R_CC_AIWCM = abs(fun_CC(Train,R_AIWCM,R_KA));
    %%DWCM
%     [R_DWCM,W] = (fun_DWCM(Train,R_KA,x0));
    %%MKA+CC=KCC
%     R_MKA = abs(fun_MKA(Train,R_KA,x0));
    R_MKA = abs(fun_CC(Train,R_NSCM,R_KA));
    R_KCC_AIWCM = abs(fun_CC(Train,R_SCM,R_MKA));
    %%%误差比较
    error_SCM(i) = sum(sum(sqrt(abs(R_SCM-rouR).^2)))/ALL;
    error_NSCM(i) = sum(sum(sqrt(abs(R_NSCM-rouR).^2)))/ALL;
    error_AML(i) = sum(sum(sqrt(abs(R_AML-rouR).^2)))/ALL;
    error_AIWCM(i) = sum(sum(sqrt(abs(R_AIWCM-rouR).^2)))/ALL;
    error_CC_SCM(i) = sum(sum(sqrt(abs(R_CC_SCM-rouR).^2)))/ALL;
    error_CC_NSCM(i) = sum(sum(sqrt(abs(R_CC_NSCM-rouR).^2)))/ALL;
    error_CC_AML(i) = sum(sum(sqrt(abs(R_CC_AML-rouR).^2)))/ALL;
    error_CC_AIWCM(i) = sum(sum(sqrt(abs(R_CC_AIWCM-rouR).^2)))/ALL;
%     error_DWCM(i) = sum(sum(sqrt(abs(R_DWCM-rouR).^2)))/ALL;
    error_KCC_AIWCM(i) = sum(sum(sqrt(abs(R_KCC_AIWCM-rouR).^2)))/ALL;
    error_MKA(i) = sum(sum(sqrt(abs(R_MKA-rouR).^2)))/ALL;
    error_KA(i) = sum(sum(sqrt(abs(R_KA-rouR).^2)))/ALL;
end
close(h)
%%误差均值
mean_error_SCM = mean(error_SCM);
mean_error_NSCM = mean(error_NSCM);
mean_error_AML = mean(error_AML);
mean_error_AIWCM = mean(error_AIWCM);
mean_error_CC_SCM = mean(error_CC_SCM);
mean_error_CC_NSCM = mean(error_CC_NSCM);
mean_error_CC_AML = mean(error_CC_AML);
mean_error_CC_AIWCM = mean(error_CC_AIWCM);
mean_error_KCC_AIWCM = mean(error_KCC_AIWCM);
mean_error_MKA = mean(error_MKA);
mean_error_KA = mean(error_KA);
% mean_error_DWCM = mean(error_DWCM);

%%误差方差
var_error_SCM = var(error_SCM);
var_error_NSCM = var(error_NSCM);
var_error_AML = var(error_AML);
var_error_AIWCM = var(error_AIWCM);
var_error_CC_SCM = var(error_CC_SCM);
var_error_CC_NSCM = var(error_CC_NSCM);
var_error_CC_AML = var(error_CC_AML);
var_error_CC_AIWCM = var(error_CC_AIWCM);
var_error_KCC_AIWCM = var(error_KCC_AIWCM);
%%画图
figure(1)
subplot(4,2,1)
plot(error_CC_NSCM,'g')
legend('CCNSCM')
axis([0,MC,0,1])
subplot(4,2,2)
plot(error_SCM,'r')
legend('SCM')
axis([0,MC,0,1])
subplot(4,2,3)
plot(error_NSCM,'b')
legend('NSCM')
axis([0,MC,0,1])
subplot(4,2,4)
plot(error_AML,'k')
legend('AML')
axis([0,MC,0,1])
subplot(4,2,5)
plot(error_CC_SCM,'y')
legend('CCSCM')
axis([0,MC,0,1])
subplot(4,2,6)
plot(error_CC_AML,'c')
legend('CCAML')
axis([0,MC,0,1])
subplot(4,2,7)
plot(error_AIWCM,'m')
legend('AIWCM')
axis([0,MC,0,1])
subplot(4,2,8)
plot(error_CC_AIWCM,'c')
legend('CCAIWCM')
axis([0,MC,0,1])
hold off


figure(2)
Xbar = {'AML','AIWCM','CCAIWCM','CCAML','CCNSCM','CCSCM','KCC','NSCM','SCM'};
Ybar = [mean_error_AML,mean_error_AIWCM,mean_error_CC_AIWCM,mean_error_CC_AML,...
        mean_error_CC_NSCM,mean_error_CC_SCM,mean_error_KCC_AIWCM,mean_error_NSCM,mean_error_SCM];
bar(Ybar,0.4)
title('均值')
set(gca,'xticklabel',Xbar)
set(gcf,'Position',[400 200 700 439])



figure(3)
Xbar = {'AML','AIWCM','CCAIWCM','CCAML','CCNSCM','CCSCM','KCC','NSCM','SCM'};
Ybar = [var_error_AML;var_error_AIWCM;var_error_CC_AIWCM;var_error_CC_AML;...
        var_error_CC_NSCM;var_error_CC_SCM;var_error_KCC_AIWCM;var_error_NSCM;var_error_SCM];
bar(Ybar,0.4)
title('方差')
set(gca,'xticklabel',Xbar)
set(gcf,'Position',[300 200 777 439])


