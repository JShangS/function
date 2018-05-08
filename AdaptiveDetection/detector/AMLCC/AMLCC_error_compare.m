clc
clear 
close all
% warning off
n = 1; %几倍的样本
str_train = 'p';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
sigma_t = [0.01,0.1:0.1:10];
% sigma_t = [11:101];
L_s = length(sigma_t);
L_R = 10000;
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
m_errorRCC = zeros(1,L_s);
m_errorRCCIter = zeros(1,L_s);
m_errorRCCML= zeros(1,L_s);
m_errorR_2 = zeros(1,L_s);
m_errorRNSCM = zeros(1,L_s);
m_errorRKA = zeros(1,L_s);
for i_s = 1:L_s
    i_s
    R_KA = zeros(size(R));
    for i = 1:1000
        t = normrnd(1,sigma_t(i_s),N,1);%%0~0.5%%失配向量
        R_KA = R_KA + R.*(t*t')/1000;
    end
    error_RCC = zeros(1,L_R);
    error_RCCIter = zeros(1,L_R);
    error_RCCML = zeros(1,L_R);
    error_R_2 = zeros(1,L_R);
    error_RNSCM = zeros(1,L_R);
    parfor i =1:L_R
        Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
        R_SCM = (fun_SCMN(Train));
        R_NSCM = (fun_NSCMN(Train));
        R_x0 = (fun_SCMN(x0));
        R_AML = fun_AML(Train);
        [R_CC,alpha(i)]=fun_CC(Train,R_SCM,R_KA);
        if  sigma_t(i_s)<0.9
            [R_CCIter,alpha_iter(i)]=fun_CCIter2(Train,R_SCM,R_KA);
        else
            [R_CCIter,alpha_iter(i)]=fun_CCIter(Train,R_SCM,R_KA);
        end
        if sigma_t(i_s)<0.9
            [R_AMLCC,alpha_aml(i)]=fun_AMLCC5(Train,R_KA);
        else
           [R_AMLCC,alpha_aml(i)]=fun_CCIter(Train,R_AML,R_KA);
        end
%         
        [R_CCML,alpha_ML(i)]=fun_MLalpha(Train,R_SCM,R_KA,x0);
        
%         R_2 = 0.5 * R_KA + 0.5 * R_NSCM;
        error_RCC(i) = norm(R_CC-R,'fro');
        error_RCCIter(i) = norm(R_CCIter-R,'fro');
        error_RAMLCC(i) = norm(R_AMLCC-R,'fro');
        error_RCCML(i) = norm(R_CCML-R,'fro');    
%         error_R_2(i) = norm(R_2-R,'fro');
        error_RNSCM(i) = norm(R_NSCM-R,'fro');
        error_RSCM(i) = norm(R_SCM-R,'fro');
        error_RAML(i) = norm(R_AML-R,'fro');
    end

    m_errorRCC(i_s) = mean(error_RCC)/norm(R,'fro')*100;
    m_errorRCCIter(i_s) = mean(error_RCCIter)/norm(R,'fro')*100;
    m_errorRAMLCC(i_s) = mean(error_RAMLCC)/norm(R,'fro')*100;
    m_errorRCCML(i_s) = mean(error_RCCML)/norm(R,'fro')*100;
%     m_errorR_2(i_s) = mean(error_R_2)/norm(R,'fro')*100;
    m_errorRNSCM(i_s) = mean(error_RNSCM)/norm(R,'fro')*100;
    m_errorRSCM(i_s) = mean(error_RSCM)/norm(R,'fro')*100;
    m_errorRAML(i_s) = mean(error_RAML)/norm(R,'fro')*100;
%     m_errorRKA(i_s) = norm(R_KA-R,'fro')/norm(R,'fro')*100;

    % mean_alpha(i_s) = mean(alpha);
    % mean_alpha_iter(i_s) = mean(alpha_iter);
    % mean_alpha_ML(i_s) = mean(alpha_ML);
end
figure
hold on 
plot(sigma_t,m_errorRNSCM,'r','LineWidth',2)%k--
plot(sigma_t,m_errorRSCM,'y','LineWidth',2)%k--
plot(sigma_t,m_errorRAML,'m','LineWidth',2)%k--
plot(sigma_t,m_errorRCC,'b','LineWidth',2)%k-.
plot(sigma_t,m_errorRAMLCC,'c','LineWidth',2)
plot(sigma_t,m_errorRCCIter,'k','LineWidth',2)
plot(sigma_t,m_errorRCCML,'g','LineWidth',2)%k:

grid on
h_leg = legend('NSCM','SCM','AML','CC','AMLCC','CCIter','CCML');
xlabel('\sigma^2','FontSize',10)
ylabel('Error','FontSize',10)
set(gca,'FontSize',10)
set(h_leg,'Location','SouthEast')
% axis([0,10,0,100])


% figure
% indexs = [1:9,20:101];%
% hold on 
% plot(sigma_t(indexs),m_errorRSCM(indexs),'k--','LineWidth',2)
% plot(sigma_t(indexs),m_errorRCC(indexs),'k-.','LineWidth',2)
% % plot(sigma_t,m_errorRCCIter_new,'b','LineWidth',2)
% plot(sigma_t(indexs),m_errorRCCIter(indexs),'k','LineWidth',2)
% plot(sigma_t(indexs),m_errorRCCML(indexs),'k:','LineWidth',2)
% % plot(sigma_t,m_errorR_2,'y','LineWidth',2)
% grid on