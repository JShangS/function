clc
clear 
close all
Read_Display_Data
Data_process
load(matFile) 
%%%%参数设置
alpha_GLC = 0;
alpha_CC=0;
beta_GLC=0;
LL=10000;
n = 3; %几倍的样本
str_train = 'p';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
sigma_t=0.1;
lambda = 1.4495;
mu = 1.3820;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
rou = 0.95;  %%协方差矩阵生成的迟滞因子
%%%Pd_CLGLRT_2Kmu1lambda3s0.1o1_p：2K：训练单元数目，mu，lambda，s：失配向量方差，
%%o1:opt=1，p：IG纹理复合高斯
%%%%假设参数设置
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
SNRout=0; % 输出SNR
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
rouR = zeros(N,N);  %%真实的杂波协方差
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
for i=1:N
    for j=1:N
        rouR(i,j)=rou^abs(i-j);%*exp(1j*2*pi*abs(i-j)*theta_sig);
    end
end
irouR=inv(rouR);
rouR_abs=abs(rouR);
% R_KA = zeros(size(rouR));
% sigma_t = 0.1;
% t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
% R_KA = R_KA+rouR.*(t*t');
% load R_KA_19980223_170435.mat
% t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
% R_KA = rouR.*(t*t');
iR_KA = inv(R_KA);
rouR_half=rouR^0.5;
Zhh = sig;
for ii = 1:LL
    ii
% Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
% x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
index_t1 = ceil(rand()*(M-10));
Train1 = Zhh(index_t1:index_t1+7,Range-L/2:Range-1);
Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2);
Train = [Train1,Train2];%%产生的训练数据,协方差矩阵为rouR的高斯杂波
x0 = Zhh(index_t1:index_t1+N-1,Range) ; % 接收信号仅包括杂波和噪声
R_SCM = (fun_SCM(Train));
R_SCMN = (fun_SCMN(Train));
iR_SCMN = inv(R_SCMN);
alpha=sqrt(SNRnum/abs(s'*irouR*s)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
x0=alpha*s+x0;
%     alpha_GLC_t = x0'*iR_KA*x0/N;
%     beta_GLC_t =  x0'*iR_SCMN*x0/N;
    [ t1,alpha_GLC_t,t2 ] = fun_CLGLRT_icovariance(R_KA,R_SCMN,x0,s,1 );
    [~,alpha_CC_t] =  fun_CC(Train,R_SCMN,R_KA);
    alpha_GLC=alpha_GLC+alpha_GLC_t/LL;
%     beta_GLC=beta_GLC_t+beta_GLC_t/LL;
    alpha_CC=alpha_CC+alpha_CC_t/LL;
end
alpha_GLC
% beta_GLC
alpha_CC