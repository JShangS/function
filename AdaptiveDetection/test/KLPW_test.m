clc
clear 
close all
% warning off
n = 2; %几倍的样本
str_train = 'g';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
sigma_t = 0.1;
opt = 'k';
rou = 0.95;  %%协方差矩阵生成的迟滞因子
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
R = fun_rho(rou,N,2);
SNRout=-5:1:20; % 输出SNR
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
R_KA = zeros(size(R));
for i = 1:1000
    t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
    R_KA = R_KA + R.*(t*t')/1000;
end
tic
parfor i =1:1000
%     i
    Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
    R_SCM = (fun_SCMN(Train));
    R_NSCM = (fun_NSCMN(Train));
    R_x0 = (fun_SCMN(x0));
    %%%%%%%%%%%%%%%%%%%%%%%%%
%     R_iter = R_x0;
%     for j = 1:iter
%         opta(i,j) = min(1,abs(trace((R_iter - R_SCM) / (R_KA - R_SCM))) / N);
%         R_iter = opta(i,j) * R_KA + (1 - opta(i,j)) * R_SCM;
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [R_result,a(i)] = fun_LogCC_new(Train,R_KA);
%     [R_result,a(i)] = fun_CLAML2(Train,R_KA);
    [R_CC,alpha(i)]=fun_CC(Train,R_SCM,R_KA);
    [R_CCIter,alpha_iter(i,:)]=fun_CCIter(Train,R_SCM,R_KA);
    [R_CCML,alpha_ML(i)]=fun_MLalpha(Train,R_SCM,R_KA,x0);
    R_AML = fun_AML(Train);
    R_2 = 0.5 * R_KA + 0.5 * R_SCM;
%     error_R(i) = norm(R_result-R,'fro');
    error_RCC(i) = norm(R_CC-R,'fro');
    error_RCCIter(i) = norm(R_CCIter-R,'fro');
    error_RCCML(i) = norm(R_CCML-R,'fro');
    error_R_2(i) = norm(R_2-R,'fro');
    error_RSCM(i) = norm(R_SCM-R,'fro');
    error_RNSCM(i) = norm(R_NSCM-R,'fro');
    error_RAML(i) = norm(R_AML-R,'fro');
end
toc
% m_errorR = mean(error_R)/norm(R,'fro')*100;
m_errorRCC = mean(error_RCC)/norm(R,'fro')*100;
m_errorRCCIter = mean(error_RCCIter)/norm(R,'fro')*100;
m_errorRCCML = mean(error_RCCML)/norm(R,'fro')*100;
m_errorR_2 = mean(error_R_2)/norm(R,'fro')*100;
m_errorRSCM = mean(error_RSCM)/norm(R,'fro')*100;
m_errorRNSCM = mean(error_RNSCM)/norm(R,'fro')*100;
m_errorAML = mean(error_RAML)/norm(R,'fro')*100;
m_errorRKA = norm(R_KA-R,'fro')/norm(R,'fro')*100;

% mean_a = mean(a);
mean_alpha = mean(alpha);
mean_alpha_iter = mean(alpha_iter);
mean_alpha_ML = mean(alpha_ML);

