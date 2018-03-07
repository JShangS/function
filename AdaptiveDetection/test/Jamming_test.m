%%%辅助数据协方差中加入含有目标信号的干扰
clc
clear 
close all
rou = 0.95;
rouj = 0.94;
clutter = 'gamma';%%杂波类型
lambda = 3; %形状参数
mu = 0.5; %尺度参数
opt = 1;
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
L = round(2*N); 
R = fun_rho(rou,N);
iR = inv(R);
PFA = 1e-3;% PFA=1e-4;
SNRout = 0:1:25; % 输出SNR
SNRnum = 10.^(SNRout/10);
MonteCarloPfa = 1/PFA*100;
MonteCarloPd = 1e4;
theta_sig = 0.3;%%系统归一化多普勒
theta_jam = 0.1;%%干扰归一化多普勒
nn = 0:N-1;
vt = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
vj = exp(-1i*2*pi*nn*theta_jam)'; %%%%%% 干扰导向矢量
%%干扰
JNR = 15;%%干燥比15dB
JNR_dB = 10.^(JNR/10);
alpha_j = sqrt(JNR_dB/abs(vj'*iR*vj)); 
jamming = alpha_j*vj;
R_j = (jamming * jamming');%干扰协方差.*fun_rho(rouj,N)
%%计算检测门限
Tace_jam = zeros(1,MonteCarloPfa); %%辅助数据含有干扰时的门限
Tace = zeros(1,MonteCarloPfa); %%辅助数据不含有目标时的门限
tic
parfor i=1:MonteCarloPfa
    Train = fun_TrainData(clutter,N,L,R,lambda, mu, opt);
    R_c = fun_SCMN(Train);
    x = fun_TrainData(clutter,N,1,R,lambda, mu, opt);
    Tace(i) = fun_ANMF(R_c,x,vt);
    Train_z = fun_TrainData(clutter,N,L,R+R_j,lambda, mu, opt);
%     R_z = R_c + R_j;
    R_z = fun_SCMN(Train_z);
    Tace_jam(i) = fun_ANMF(R_z,x,vt);
    Topt(i) = abs(vt'*iR*x)^2/(vt'*iR*vt)/(x'*iR*x);
end
toc
TACE=sort(Tace,'descend');
TACE_JAM=sort(Tace_jam,'descend');
TOPT=sort(Topt,'descend');
Th_OPT=(TOPT(floor(MonteCarloPfa*PFA-1))+TOPT(floor(MonteCarloPfa*PFA)))/2;
Th_ACE=(TACE(floor(MonteCarloPfa*PFA-1))+TACE(floor(MonteCarloPfa*PFA)))/2;
Th_ACE_jam=(TACE_JAM(floor(MonteCarloPfa*PFA-1))+TACE_JAM(floor(MonteCarloPfa*PFA)))/2;
%%计算检测概率
counter_ace = 0;
counter_opt = 0; 
counter_ace_jam = 0;
counter_ace_orig = 0;

Pd_ACE_mc = zeros(1,length(SNRout));
Pd_ACE_jam_mc = zeros(1,length(SNRout));
Pd_ACE_orig_mc = zeros(1,length(SNRout));
Pd_OPT_mc = zeros(1,length(SNRout));
alpha=sqrt(SNRnum/abs(vt'*iR*vt)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
h = waitbar(0,'Please wait...');
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
    parfor i=1:MonteCarloPd
        Train = fun_TrainData(clutter,N,L,R,lambda, mu, opt);
        R_c = fun_SCMN(Train);
        Train_z = fun_TrainData(clutter,N,L,R+R_j,lambda, mu, opt);
%         R_z = R_c + R_j;
        R_z = fun_SCMN(Train_z);
        x = fun_TrainData(clutter,N,1,R,lambda, mu, opt);
        x = alpha(m)*vt+x;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
        Tace = fun_ANMF(R_z,x,vt); %%辅助数据协方差中含有干扰的ACE
        Tace_orig = fun_ANMF(R_c,x,vt); %%原始的ACE
        Topt = abs(vt'*iR*x)^2/(vt'*iR*vt)/(x'*iR*x);
        if Tace > Th_ACE;         counter_ace = counter_ace + 1;          end
        if Tace > Th_ACE_jam;         counter_ace_jam = counter_ace_jam + 1;          end
        if Tace_orig > Th_ACE;    counter_ace_orig = counter_ace_orig + 1;          end
        if Topt > Th_OPT;          counter_opt = counter_opt+1;           end
    end
    Pd_ACE_mc(m) = counter_ace/MonteCarloPd;          counter_ace = 0;
    Pd_ACE_jam_mc(m) = counter_ace_jam / MonteCarloPd;   counter_ace_jam = 0; 
    Pd_ACE_orig_mc(m) = counter_ace_orig / MonteCarloPd;   counter_ace_orig = 0; 
    Pd_OPT_mc(m) = counter_opt / MonteCarloPd;   counter_opt = 0; 
end
close(h)
figure()
hold on
% plot(SNRout,Pd_ACE_mc,'r');
plot(SNRout,Pd_ACE_jam_mc,'k');
plot(SNRout,Pd_ACE_orig_mc,'b');
plot(SNRout,Pd_OPT_mc,'r');
legend('(门限确定也含干扰)含干扰的ACE','没干扰的ACE','最优')
% legend('含干扰的ACE','门限确定含有干扰的ACE','没干扰的ACE','最优')
str = ['信号归一化多普勒 = ', num2str(theta_sig),',干扰归一化多普勒 = ', num2str(theta_jam)];
title(str);