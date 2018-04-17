clc
clear 
close all
% warning off
n = 0.5; %����������
str_train = 'p';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
opt = 'k';
rou = 0.95;  %%Э����������ɵĳ�������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
R = fun_rho(rou,N,1,0.2);
SNRout=-5:1:20; % ���SNR
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% ϵͳ����ʸ��
iter = 10;
tic
for i =1:100
    i
    Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
    R_SCM = (fun_SCMN(Train));
    R_NSCM = (fun_NSCMN(Train));
    R_x0 = (fun_SCMN(x0));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    R_LogE = (fun_RLogEMean(Train));
    R_LogEP = (fun_RLogMeanPer(Train));
    R_P = fun_Perm(Train);

    error_RLogE(i) = norm(R_LogE - R,'fro');
    error_RLogEP(i) = norm(R_LogEP - R,'fro');
    error_RP(i) = norm(R_P - R,'fro');
    error_RSCM(i) = norm(R_SCM - R,'fro');
end
toc
m_errorRLogE = mean(error_RLogE);
m_errorRLogEP = mean(error_RLogEP);
m_errorRP = mean(error_RP);
m_errorRSCM = mean(error_RSCM);

