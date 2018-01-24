clc
clear 
close all
%%%%��������
alpha_GLC = 0;
alpha_CC=0;
beta_GLC=0;
LL=10000;
n = 1; %����������
str_train = 'p';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
rou = 0.95;  %%Э����������ɵĳ�������
%%%Pd_CLGLRT_2Kmu1lambda3s0.1o1_p��2K��ѵ����Ԫ��Ŀ��mu��lambda��s��ʧ���������
%%o1:opt=1��p��IG�����ϸ�˹
%%%%�����������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
SNRout=0; % ���SNR
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
rouR = zeros(N,N);  %%��ʵ���Ӳ�Э����
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% ϵͳ����ʸ��
for i=1:N
    for j=1:N
        rouR(i,j)=rou^abs(i-j);%*exp(1j*2*pi*abs(i-j)*theta_sig);
    end
end
irouR=inv(rouR);
rouR_abs=abs(rouR);
R_KA = zeros(size(rouR));
sigma_t = 0.5;
t = normrnd(1,sigma_t,N,1);%%0~0.5%%ʧ������
R_KA = R_KA+rouR.*(t*t');
% load R_KA_0.1.mat
iR_KA = inv(R_KA);
rouR_half=rouR^0.5;
for ii = 1:LL
    ii
Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % �����źŽ������Ӳ�������
R_SCM = (fun_SCM(Train));
R_SCMN = (fun_SCMN(Train));
iR_SCMN = inv(R_SCMN);
alpha=sqrt(SNRnum/abs(s'*irouR*s)); % ����SNR=|alpha|^2*s'*R^(-1)*s���|alpha|
x0=alpha*s+x0;
%     alpha_GLC_t = x0'*iR_KA*x0/N;
%     beta_GLC_t =  x0'*iR_SCMN*x0/N;
    [ t1,alpha_GLC_t,t2 ] = fun_CLGLRT_covariance2(R_KA,R_SCMN,x0,s,1 );
    [~,alpha_CC_t] =  fun_CC(Train,R_SCMN,R_KA);
    alpha_GLC=alpha_GLC+alpha_GLC_t/LL;
%     beta_GLC=beta_GLC_t+beta_GLC_t/LL;
    alpha_CC=alpha_CC+alpha_CC_t/LL;
end
alpha_GLC
% beta_GLC
alpha_CC