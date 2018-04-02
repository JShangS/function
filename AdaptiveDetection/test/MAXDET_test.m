clc
clear 
close all
echo on
rou = 0.95;
rouM = [0.1,0.5,0.9];
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
L = round(2*N); 
R = fun_rho(rou,N);
RM =  fun_rho(rouM,N);
iR = inv(R);
t = normrnd(1,0.1,N,1);%%0~0.5%%ʧ������
M0 = inv(R.*(t*t'));
PFA = 1e-3;% PFA=1e-4;
SNRout = 0:1:20; % ���SNR
SNRnum = 10.^(SNRout/10);
MonteCarloPfa = 1/PFA*100;
MonteCarloPd = 1e4;
theta_sig = 0.1;%%ϵͳ��һ��������
theta_jam = 0.15;%%���Ź�һ��������
nn = 0:N-1;
vt = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% ϵͳ����ʸ��
Train = fun_TrainData('g',N,L,R,3);
S = abs(fun_NSCMN(Train));
R_DET = sdpvar(N,N);
t = sdpvar(1,3);
F = set( R_DET > 0); %��������
F = F + set(t(1) * RM(:,:,1) + t(2) * RM(:,:,2) + t(3) * RM(:,:,3) - R_DET == 0);
solution = solvesdp(F,-logdet(R_DET)+trace(R_DET*S));

R_DET = inv(double(R_DET));
t = double(t);
