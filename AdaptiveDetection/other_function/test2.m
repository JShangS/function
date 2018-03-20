clc
clear 
close all
% rou = 1;  %%Э����������ɵĳ�������
% N = 8;
% theta_sig = 0.1;
% for i=1:N
%     for j=1:N
%         rouR(i,j)=rou^(i-j)*exp(1j*2*pi*(i-j)*theta_sig);%
%     end
% end
n = 4; %����������
str_train = 'g';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
opt_train = 3; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
sigma_t = 1;
rou = 0.95;  %%Э����������ɵĳ�������
%%%Pd_CLGLRT_2Kmu1lambda3s0.1o1_p��2K��ѵ����Ԫ��Ŀ��mu��lambda��s��ʧ���������
%%o1:opt=1��p��IG�����ϸ�˹
%%%%�����������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% ϵͳ����ʸ��
L=round(n*N); 
for i=1:N
    for j=1:N
        rouR(i,j)=rou^abs(i-j);%*exp(1j*2*pi*abs(i-j)*theta_sig);
    end
end
Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
R_SCM = (fun_SCMN(Train));
R_NSCM = fun_NSCMN(Train);
R_AML = fun_AML(Train);
M = 1000;
ft = linspace(-0.5,0.5,M);
PSD=fun_PSD(rouR,ft);
% PSD=PSD/max(abs(PSD));
% PSD=fun_value2dB(PSD);
plot(ft,(PSD))
hold on
plot(ft,abs(fun_PSD(R_SCM,ft)),'r')
plot(ft,(fun_PSD(R_NSCM,ft)),'k.-')
plot(ft,(fun_PSD(R_AML,ft)),'g.-')
