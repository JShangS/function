clc
clear 
close all
echo on
rou = 0.95;
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
L = round(2*N); 
R = fun_rho(rou,N);
iR = inv(R);
t = normrnd(1,0.1,N,1);%%0~0.5%%失配向量
M0 = inv(R.*(t*t'));
PFA = 1e-3;% PFA=1e-4;
SNRout = 0:1:20; % 输出SNR
SNRnum = 10.^(SNRout/10);
MonteCarloPfa = 1/PFA*100;
MonteCarloPd = 1e4;
theta_sig = 0.1;%%系统归一化多普勒
theta_jam = 0.15;%%干扰归一化多普勒
nn = 0:N-1;
vt = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
Train = fun_TrainData('g',N,L,R,3);
S = fun_SCMN(Train);
R_DET = sdpvar(N,N);
F = set( R_DET > 0); %限制条件
solution = solvesdp(F,-logdet((R_DET))+trace(R_DET*S));

