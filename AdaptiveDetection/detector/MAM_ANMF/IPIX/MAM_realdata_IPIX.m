%%%实现一个基于色加载的GLRT检测器
%%%%
%%%%多线程的时候不要在跑程序的时候改程序!!!!!
clc
clear 
close all
Read_Display_Data
Data_process
load(matFile) 
%%%%参数设置
n = 2; %几倍的样本
sigma_t = 0.0;
% range = 14;
%%参数估计来自《Maximum Likelihood Estimation for
%%%            Compound-Gaussian Clutter with Inverse GammaTexture》

lambda = 3;
mu = 1;
% lambda = 1.4495;
% mu = 1.3820;
str_train = 'g';
opt_train = 1;
%%%%参数设置
% N = 8;
SNRout=-5:1:25; % 输出SNR
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=100/PFA;
MonteCarloPd=1e4;
rouR = R_KA;  %%真实的杂波协方差
irouR = inv(rouR);
rouM=[0.9,0.95,0.99];%%%%%%%%%MAM模型
for i =1:length(rouM)
    R = fun_rho(rouM(i),N);
    MAM(:,:,i)=R;
end
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
%%%%%正式开始%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%门限计算%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = waitbar(0,'Please wait...');
Zhh = sig;
tic
parfor i = 1:MonteCarloPfa
%     waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
%%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
%%%%%%%%%%%%%用杂波模型确定门限%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
%       x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
%%%%%%%%用实测数据确定门限，但是数据有限。%%%%%%%%%%%%%%%%%%%%%%%
    index_t1 = ceil(rand()*(M-10));
    Train1 = Zhh(index_t1:index_t1+7,Range-L/2:Range-1);
    Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2);
    Train = [Train1,Train2];%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = Zhh(index_t1:index_t1+N-1,Range) ; % 接收信号仅包括杂波和噪声
    %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%
    R_SCM = (fun_SCM(Train));
    
    R_SCMN = (fun_SCMN(Train));
    
    R_NSCM = (fun_NSCM(Train));
    
    R_NSCMN = (fun_NSCMN(Train));
    %%%%%%%MAM协方差估计%%%%%%%%%%%%%
    Rx0=fun_SCM(x0);
    R_MAM = fun_information_estimation(Rx0,MAM);
    %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% ANMF_SCM
    Tanmf_SCM(i) = fun_ANMF(R_SCM,x0,s);
    %%%%%% ANMF_NSCM
    Tanmf_NSCM(i) = fun_ANMF(R_NSCM,x0,s);
    %%%%%% ANMF_MAM
    Tanmf_MAM(i) = fun_ANMF(R_MAM,x0,s);
    %%%%%% GLRT_MAM
    Tglrt_mam(i) = fun_MAM_GLRT(R_MAM,x0,s,lambda,mu);
end
toc
% close(h)
TANMF_SCM=sort(Tanmf_SCM,'descend');
TANMF_NSCM=sort(Tanmf_NSCM,'descend');
TANMF_MAM=sort(Tanmf_MAM,'descend');
TGLRT_MAM=sort(Tglrt_mam,'descend');

Th_SCM=(TANMF_SCM(floor(MonteCarloPfa*PFA-1))+TANMF_SCM(floor(MonteCarloPfa*PFA)))/2;
Th_NSCM=(TANMF_NSCM(floor(MonteCarloPfa*PFA-1))+TANMF_NSCM(floor(MonteCarloPfa*PFA)))/2;
Th_MAM=(TANMF_MAM(floor(MonteCarloPfa*PFA-1))+TANMF_MAM(floor(MonteCarloPfa*PFA)))/2;
Th_GLRTMAM=(TGLRT_MAM(floor(MonteCarloPfa*PFA-1))+TGLRT_MAM(floor(MonteCarloPfa*PFA)))/2;
%%%%%%%%%%%%%%%%%%%%%检测概率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

counter_scm=0;
counter_nscm=0;
counter_mam=0;
counter_glrtmam=0;

Pd_SCM_mc = zeros(1,length(SNRout));
Pd_NSCM_mc = zeros(1,length(SNRout));
Pd_MAM_mc = zeros(1,length(SNRout));
Pd_GLRTMAM_mc = zeros(1,length(SNRout));
alpha=sqrt(SNRnum/abs(s'*irouR*s)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
% alpha=sqrt(SNRnum/abs(s'*irouR*s));
h = waitbar(0,'Please wait...');
tic
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
    parfor i=1:MonteCarloPd 
 %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
%%%%%%%%%%%%%用杂波模型产生训练数据%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
%         x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        index_t1 = ceil(rand()*(M-10));
        Train1 = Zhh(index_t1:index_t1+7,Range-L/2:Range-1);
        Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2);
        Train = [Train1,Train2];%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        x0 = Zhh(index_t1:index_t1+7,Range) ; % 接收信号仅包括杂波和噪声
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x0=alpha(m)*s+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
        %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%
        R_SCM = (fun_SCM(Train));
    
        R_SCMN = (fun_SCMN(Train));
    
        R_NSCM = (fun_NSCM(Train));
    
        R_NSCMN = (fun_NSCMN(Train));
        %%%%%%%MAM协方差估计%%%%%%%%%%%%%
        Rx0=fun_SCM(x0);
        R_MAM = fun_information_estimation(Rx0,MAM);     
        %%%检测信号
        x0=alpha(m)*s+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
        %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% ANMF_SCM
        Tscm = fun_ANMF(R_SCM,x0,s);
        %%%%%% ANMF_NSCM
        Tnscm = fun_ANMF(R_NSCM,x0,s);
        %%%%%% ANMF_MAM
        Tmam = fun_ANMF(R_MAM,x0,s);
        %%%%%% GLRT_MAM
        Tglrtmam = fun_MAM_GLRT(R_MAM,x0,s,lambda,mu);
        %%%判断%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        if Tscm>Th_SCM;          counter_scm=counter_scm+1;        end                
        if Tnscm>Th_NSCM;       counter_nscm=counter_nscm+1;    end   
        if Tmam>Th_MAM;      counter_mam=counter_mam+1;    end
        if Tglrtmam>Th_GLRTMAM;      counter_glrtmam=counter_glrtmam+1;    end
    end
    Pd_SCM_mc(m)=counter_scm/MonteCarloPd;           counter_scm=0;
    Pd_NSCM_mc(m)=counter_nscm/MonteCarloPd;        counter_nscm=0;
    Pd_MAM_mc(m)=counter_mam/MonteCarloPd;       counter_mam=0; 
    Pd_GLRTMAM_mc(m)=counter_glrtmam/MonteCarloPd;       counter_glrtmam=0;
end
close(h)
toc
figure(1);
hold on
plot(SNRout,Pd_SCM_mc,'b-+','linewidth',2)
plot(SNRout,Pd_NSCM_mc,'k.-','linewidth',2)
plot(SNRout,Pd_MAM_mc,'g-s','linewidth',2)
plot(SNRout,Pd_GLRTMAM_mc,'y-o','linewidth',2)
h_leg = legend('ANMF with SCM','ANMF with NSCM','ANMF with MAM','GLRT with MAM');
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
grid on
% str=['Pd_MAM_IPIX_',num2str(n),'K','_',cdfFile_t,'.mat'];
% save(str,'SNRout','Pd_SCM_mc','Pd_NSCM_mc',...
%          'Pd_MAM_mc');