clc
clear 
close all
% warning off
n = 2; %����������
str_train = 'p';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
sigma_t = 0.1;
rou = 0.95;  %%Э����������ɵĳ�������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
R = fun_rho(rou,N,2);
SNRout=-5:1:20; % ���SNR
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% ϵͳ����ʸ��
R_KA = zeros(size(R));
for i = 1:1000
    t = normrnd(1,sigma_t,N,1);%%0~0.5%%ʧ������
    R_KA = R_KA + R.*(t*t')/1000;
end
iter = 10;
for i =1:1000
    i
    Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
    R_SCM = abs(fun_SCMN(Train));
    R_x0 = abs(fun_SCMN(x0));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    [R_result,a(i)] = fun_LogEDCC(Train,R_KA);
    [R_CC,alpha(i)]=fun_CC(Train,R_SCM,R_KA);
    R_2 = 0.5 * R_KA + 0.5 * R_SCM;
    error_R(i) = norm(R_result-R,'fro')/norm(R,'fro');
    error_RCC(i) = norm(R_CC-R,'fro')/norm(R,'fro');
    error_R_2(i) = norm(R_2-R,'fro')/norm(R,'fro');
    error_RSCM(i) = norm(R_SCM-R,'fro')/norm(R,'fro');
end
m_errorR = mean(error_R)*100;
m_errorRCC = mean(error_RCC)*100;
m_errorR_2 = mean(error_R_2)*100;
m_errorRSCM = mean(error_RSCM)*100;
m_errorRKA = norm(R_KA-R,'fro')/norm(R,'fro')*100;
mean_a = mean(a);
mean_alpha = mean(alpha);
plot(a,'r')
hold on
plot(alpha,'k')
axis([1,1000,0,1]);

