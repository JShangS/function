%%%����Э������ƵıȽ�
%%%���������õ�norm��A,'fro'��
clc
clear 
close all
Na = 4;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
rou = 0.95;  %%Э����������ɵĳ�������
rouR = zeros(N,N);  %%��ʵ���Ӳ�Э����
L=round(2*N); 
% theta_sig = 0.1;
% nn = 0:N-1;
% s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% ϵͳ����ʸ��
for i=1:N
    for j=1:N
        rouR(i,j)=rou^abs(i-j);
    end
end
irouR=inv(rouR);
sigma_t = 0.3;
t = normrnd(1,sigma_t,N,1);%%0~0.5%%ʧ������
R_KA = rouR.*(t*t');
rouR_half=rouR^0.5;
MC = 1000;
ALL = norm(rouR,'fro');%����Э�����Frobenius����
error_SCM = zeros(MC,1);
error_NSCM = zeros(MC,1);
error_AML = zeros(MC,1);
error_AIWCM = zeros(MC,1);
error_CC_SCM = zeros(MC,1);
error_CC_NSCM = zeros(MC,1);
error_CC_AML = zeros(MC,1);
error_CC_AIWCM = zeros(MC,1);
error_KCC_AIWCM = zeros(MC,1);
error_MKA = zeros(MC,1);
error_KA = zeros(MC,1);
h = waitbar(0,'Please wait...');
for i = 1:MC
    waitbar(i/MC,h,sprintf([num2str(i/MC*100),'%%']));
    %%�����Ӳ�������
    Train = fun_TrainData('g',N,L,rouR);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    x0 = fun_TrainData('g',N,1,rouR); % �����źŽ������Ӳ�������
    Sx0 = (x0*x0');
%%Э�������
    %%SCM������Э����
    R_SCM = abs(fun_SCM(Train));
    %%NSCM����һ������Э����
    R_NSCM = abs(fun_NSCM(Train));
    %%AML�����������Ȼ
    R_AML = abs(fun_AML(Train));
    %%AIWCM������Ӧ����Э�������
    R_AIWCM =abs(fun_AIWCM(Train,x0));
    %%CC+R_SCM��͹�Ż���֪ʶ����+Э��������㷨��ѵ������Э����
    [R_CC_SCM,alpha0(i)] = (fun_CC(Train,R_SCM,R_KA));
    R_CC_SCM = abs(R_CC_SCM);
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
    R_KCC_AIWCM = abs(fun_CC(Train,R_AML,R_MKA));
    %%MLalpha%%����,����Ҫ
%     [R_MLalpha,alpha(i)] = fun_MLalpha(Train,R_SCM,R_KA,x0);
%     R_MLalpha = abs(R_MLalpha);
    %%%���Ƚ�
    error_SCM(i) = norm(R_SCM-rouR,'fro')/ALL;
    error_NSCM(i) = norm(R_NSCM-rouR,'fro')/ALL;
    error_AML(i) = norm(R_AML-rouR,'fro')/ALL;
    error_AIWCM(i) = norm(R_AIWCM-rouR,'fro')/ALL;
    error_CC_SCM(i) = norm(R_CC_SCM-rouR,'fro')/ALL;
    error_CC_NSCM(i) = norm(R_CC_NSCM-rouR,'fro')/ALL;
    error_CC_AML(i) = norm(R_CC_AML-rouR,'fro')/ALL;
    error_CC_AIWCM(i) = norm(R_CC_AIWCM-rouR,'fro')/ALL;
%     error_DWCM(i) = norm(R_DWCM-rouR,'fro')/ALL;
    error_KCC_AIWCM(i) = norm(R_KCC_AIWCM-rouR,'fro')/ALL;
    error_MKA(i) = norm(R_MKA-rouR,'fro')/ALL;
    error_KA(i) = norm(R_KA-rouR,'fro')/ALL;
%     error_MLalpha(i) = norm(R_MLalpha-rouR,'fro')/ALL;
end
close(h)
%%����ֵ
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
% mean_error_MLalpha = mean(error_MLalpha);
%%����
var_error_SCM = var(error_SCM);
var_error_NSCM = var(error_NSCM);
var_error_AML = var(error_AML);
var_error_AIWCM = var(error_AIWCM);
var_error_CC_SCM = var(error_CC_SCM);
var_error_CC_NSCM = var(error_CC_NSCM);
var_error_CC_AML = var(error_CC_AML);
var_error_CC_AIWCM = var(error_CC_AIWCM);
var_error_KCC_AIWCM = var(error_KCC_AIWCM);
% var_error_DWCM = var(error_DWCM);
% var_error_MLalpha = mean(error_MLalpha);
%%��ͼ
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
Xbar = {'AML','AIWCM','CCAIWCM','CCAML','CCNSCM','CCSCM','KCC','NSCM','SCM'};%,'ML'
Ybar = [mean_error_AML,mean_error_AIWCM,mean_error_CC_AIWCM,mean_error_CC_AML,...
        mean_error_CC_NSCM,mean_error_CC_SCM,mean_error_KCC_AIWCM,mean_error_NSCM,...
        mean_error_SCM];%mean_error_MLalpha
bar(Ybar,0.4)
title('��ֵ')
set(gca,'xticklabel',Xbar)
set(gcf,'Position',[400 200 777 439])



figure(3)
Xbar = {'AML','AIWCM','CCAIWCM','CCAML','CCNSCM','CCSCM','KCC','NSCM','SCM'};%,'ML'
Ybar = [var_error_AML;var_error_AIWCM;var_error_CC_AIWCM;var_error_CC_AML;...
        var_error_CC_NSCM;var_error_CC_SCM;var_error_KCC_AIWCM;var_error_NSCM;...
        var_error_SCM];%var_error_MLalpha
bar(Ybar,0.4)
title('����')
set(gca,'xticklabel',Xbar)
set(gcf,'Position',[300 200 777 439])


