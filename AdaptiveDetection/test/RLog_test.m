clc
clear 
close all
% warning off
n = 2; %几倍的样本
str_train = 'p';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
sigma_t = 10;
rou = 0.95;  %%协方差矩阵生成的迟滞因子
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
R = fun_rho(rou,N,1,0);
SNRout=-5:1:20; % 输出SNR
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
iter = 10;
R_KA = zeros(size(R));
for i = 1:1000
    t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
    R_KA = R_KA + R.*(t*t')/1000;
end
tic
% mu = zeros(N,1);
for i =1:100
    i
    Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
    R_SCM = (fun_SCMN(Train));
    R_NSCM = (fun_NSCMN(Train));
    R_x0 = (fun_SCMN(x0));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    R_LogE = fun_RLogEMean(Train);
    [R_LogNormCC,alpha_log(i)] = fun_LogEDCC(Train,R_NSCM,R_KA);
%     R_LogAML = fun_LogAML(Train);
    R_AML = fun_AML(Train);
%     R_LogEP = (fun_RLogMeanPer(Train));
%     R_P = fun_Perm(Train);

    error_RLogE(i) = norm(R_LogE - R,'fro');
    error_RLogNormCC(i) = norm(R_LogNormCC - R,'fro');
%     error_RLogAML(i) = norm(R_LogAML - R,'fro');
    error_RAML(i) = norm(R_AML - R,'fro');
%     error_RLogEP(i) = norm(R_LogEP - R,'fro');
%     error_RP(i) = norm(R_P - R,'fro');
    error_RSCM(i) = norm(R_SCM - R,'fro');
    error_RNSCM(i) = norm(R_NSCM - R,'fro');
end
toc
m_errorRLogE = mean(error_RLogE)/norm(R,'fro');
m_errorRLogNormCC = mean(error_RLogNormCC)/norm(R,'fro');
% m_errorRLogAML = mean(error_RLogAML)/norm(R,'fro');
m_errorRAML = mean(error_RAML)/norm(R,'fro');
% m_errorRLogEP = mean(error_RLogEP);
% m_errorRP = mean(error_RP);
m_errorRSCM = mean(error_RSCM)/norm(R,'fro');
m_errorRNSCM = mean(error_RNSCM)/norm(R,'fro');
m_alpha_log = mean(alpha_log);

