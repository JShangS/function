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
% range = 14;
%%参数估计来自《Maximum Likelihood Estimation for
%%%            Compound-Gaussian Clutter with Inverse GammaTexture》

lambda = 3;
mu = 1;
% lambda = 1.4495;
% mu = 1.3820;
% str_train = 'g';
opt_train = 1;
%%%%参数设置
% N = 8;
SNRout=-5:1:25; % 输出SNR
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPd=1e4;
MonteCarloPfa=M-MonteCarloPd;
rouR = R_KA;  %%真实的杂波协方差
irouR = inv(rouR);
rouM=[0.1,0.2,0.55,0.6,0.95,0.8];%%%%%%%%%MAM模型
for i =1:length(rouM)
    R = fun_rho(rouM(i),N);%-0.23
    MAM(:,:,i)=R;
end
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'/sqrt(N); %%%%%% 系统导向矢量
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
    index_t1 = i;
    Train1 = Zhh(index_t1:index_t1+7,Range-L/2:Range-1);
    Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2);
    Train = [Train1,Train2];%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = Zhh(index_t1:index_t1+N-1,Range) ; % 接收信号仅包括杂波和噪声
   %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%
    R_SCM = (fun_SCM(Train));
    
    R_SCMN = (fun_SCMN(Train));
    
    R_NSCM = (fun_NSCM(Train));
    
    R_NSCMN = (fun_NSCMN(Train));
    
    Rx0 = fun_Positive(x0,3);
    %%%%%%%MAM协方差估计%%%%%%%%%%%%%
%     Rx0 = x0*x0';
%     [V,D] = svd(Rx0);
%     Evalue = abs(diag(D));
%     index_1 = find(Evalue<=1);
%     Evalue(index_1) = 1;
%     % Evalue = sort(Evalue,'descend');
%     D = diag(Evalue);
%     Rx0 = abs(V)*D*abs(V');
%     if abs(det(Rx0))<1e-6 %%%%%%%%%%%处理边界值
%         continue;
%     end
%     R_MAM = MAM;%fun_congnition(Train);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R_MAM_r = fun_information_estimation(Rx0,MAM,'r');%ReimanDistance
%     R_MAM_c = fun_information_estimation(R_NSCMN,MAM,'c');%CholeskyDistance
    R_MAM_e = fun_information_estimation(Rx0,MAM,'e');%Euclidean Distance
    R_MAM_l = fun_information_estimation(Rx0,MAM,'l');%Log Euclidean Distance
    R_MAM_p = fun_information_estimation(Rx0,MAM,'p');%Power-Euclidean distance
    R_MAM_ro = fun_information_estimation(Rx0,MAM,'ro');%Root Euclidean Distance
    %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% ANMF_SCM
    Tanmf_SCM(i) = fun_AMF(R_SCMN,x0,s);
    %%%%%% ANMF_NSCM
    Tanmf_NSCM(i) = fun_ANMF(R_SCMN,x0,s);
    %%%%%% ANMF_MAM_r
    Tanmf_MAM_r(i) = fun_ANMF(R_MAM_r,x0,s);
    %%%%%% ANMF_MAM_c
%     Tanmf_MAM_c(i) = fun_ANMF(R_MAM_c,x0,s);
    %%%%%% ANMF_MAM_e
    Tanmf_MAM_e(i) = fun_ANMF(R_MAM_e,x0,s);
    %%%%%% ANMF_MAM_l
    Tanmf_MAM_l(i) = fun_ANMF(R_MAM_l,x0,s);
    %%%%%% ANMF_MAM_p
    Tanmf_MAM_p(i) = fun_ANMF(R_MAM_p,x0,s);
    %%%%%% ANMF_MAM_ro
    Tanmf_MAM_ro(i) = fun_ANMF(R_MAM_ro,x0,s);
    %%%%%% GLRT_MAM
    Tglrt_mam(i) = fun_MAM_GLRT(MAM,x0,s);
end
toc
% close(h)
TANMF_SCM=sort(Tanmf_SCM,'descend');
TANMF_NSCM=sort(Tanmf_NSCM,'descend');
TANMF_MAM_r=sort(Tanmf_MAM_r,'descend');
% TANMF_MAM_c=sort(Tanmf_MAM_c,'descend');
TANMF_MAM_e=sort(Tanmf_MAM_e,'descend');
TANMF_MAM_l=sort(Tanmf_MAM_l,'descend');
TANMF_MAM_p=sort(Tanmf_MAM_p,'descend');
TANMF_MAM_ro=sort(Tanmf_MAM_ro,'descend');
TGLRT_MAM=sort(Tglrt_mam,'descend');
MonteCarloPfa = length(TANMF_SCM);
Th_SCM=(TANMF_SCM(floor(MonteCarloPfa*PFA-1))+TANMF_SCM(floor(MonteCarloPfa*PFA)))/2;
Th_NSCM=(TANMF_NSCM(floor(MonteCarloPfa*PFA-1))+TANMF_NSCM(floor(MonteCarloPfa*PFA)))/2;

Th_MAM_r=(TANMF_MAM_r(floor(MonteCarloPfa*PFA-1))+TANMF_MAM_r(floor(MonteCarloPfa*PFA)))/2;
% Th_MAM_c=(TANMF_MAM_c(floor(MonteCarloPfa*PFA-1))+TANMF_MAM_c(floor(MonteCarloPfa*PFA)))/2;
Th_MAM_e=(TANMF_MAM_e(floor(MonteCarloPfa*PFA-1))+TANMF_MAM_e(floor(MonteCarloPfa*PFA)))/2;
Th_MAM_l=(TANMF_MAM_l(floor(MonteCarloPfa*PFA-1))+TANMF_MAM_l(floor(MonteCarloPfa*PFA)))/2;
Th_MAM_p=(TANMF_MAM_p(floor(MonteCarloPfa*PFA-1))+TANMF_MAM_p(floor(MonteCarloPfa*PFA)))/2;
Th_MAM_ro=(TANMF_MAM_ro(floor(MonteCarloPfa*PFA-1))+TANMF_MAM_ro(floor(MonteCarloPfa*PFA)))/2;

Th_GLRTMAM=(TGLRT_MAM(floor(MonteCarloPfa*PFA-1))+TGLRT_MAM(floor(MonteCarloPfa*PFA)))/2;
%%%%%%%%%%%%%%%%%%%%%检测概率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

counter_scm=0;
counter_nscm=0;
counter_mam_r=0;
% counter_mam_c=0;
counter_mam_e=0;
counter_mam_l=0;
counter_mam_p=0;
counter_mam_ro=0;
counter_glrtmam=0;

Pd_SCM_mc = zeros(1,length(SNRout));
Pd_NSCM_mc = zeros(1,length(SNRout));
Pd_MAM_r_mc = zeros(1,length(SNRout));
% Pd_MAM_c_mc = zeros(1,length(SNRout));
Pd_MAM_e_mc = zeros(1,length(SNRout));
Pd_MAM_l_mc = zeros(1,length(SNRout));
Pd_MAM_p_mc = zeros(1,length(SNRout));
Pd_MAM_ro_mc = zeros(1,length(SNRout));
Pd_GLRTMAM_mc = zeros(1,length(SNRout));
alpha=sqrt(SNRnum/abs(s'*irouR*s)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
% alpha=sqrt(SNRnum/abs(s'*irouR*s));
h = waitbar(0,'Please wait...');
tic
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
    count_Pd = MonteCarloPd;
    parfor i=1:MonteCarloPd 
 %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
%%%%%%%%%%%%%用杂波模型产生训练数据%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
%         x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        index_t1 = i+49990;
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
        
        Rx0 = fun_Positive(x0,3);
        %%%%%%%MAM协方差估计%%%%%%%%%%%%%
%         Rx0 = x0*x0';
%         [V,D] = svd(Rx0);
%         Evalue = abs(diag(D));
%         index_1 = find(Evalue<=1);
%         Evalue(index_1) = 1;
%         % Evalue = sort(Evalue,'descend');
%         D = diag(Evalue);
%         Rx0 = abs(V)*D*abs(V');
%         if abs(det(Rx0))<1e-6%%%%%%%%%%%处理边界值
%             count_Pd = count_Pd-1;
%             continue;
%         end
%         R_MAM = MAM;%fun_congnition(Train);
        R_MAM_r = fun_information_estimation(Rx0,MAM,'r');
%         R_MAM_c = fun_information_estimation(abs(Rx0),MAM,'c');
        R_MAM_e = fun_information_estimation(Rx0,MAM,'e');
        R_MAM_l = fun_information_estimation(Rx0,MAM,'l');
        R_MAM_p = fun_information_estimation(Rx0,MAM,'p');
        R_MAM_ro = fun_information_estimation(Rx0,MAM,'ro');
        %%%检测信号
        x0=alpha(m)*s+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
        %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% ANMF_SCM
        Tscm = fun_AMF(R_SCMN,x0,s);
        %%%%%% ANMF_NSCM
        Tnscm = fun_ANMF(R_SCMN,x0,s);
        %%%%%% ANMF_MAM_r
        Tmam_r = fun_ANMF(R_MAM_r,x0,s);
        %%%%%% ANMF_MAM_c
%         Tmam_c = fun_ANMF(R_MAM_c,x0,s);
        %%%%%% ANMF_MAM_e
        Tmam_e = fun_ANMF(R_MAM_e,x0,s);
        %%%%%% ANMF_MAM_l
        Tmam_l = fun_ANMF(R_MAM_l,x0,s);
        %%%%%% ANMF_MAM_p
        Tmam_p = fun_ANMF(R_MAM_p,x0,s);
        %%%%%% ANMF_MAM_ro
        Tmam_ro = fun_ANMF(R_MAM_ro,x0,s);
        %%%%%% GLRT_MAM
        Tglrtmam = fun_MAM_GLRT(MAM,x0,s);
        %%%判断%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        if Tscm>Th_SCM;          counter_scm=counter_scm+1;        end                
        if Tnscm>Th_NSCM;       counter_nscm=counter_nscm+1;    end   
        if Tmam_r>Th_MAM_r;      counter_mam_r=counter_mam_r+1;    end
%         if Tmam_c>Th_MAM_c;      counter_mam_c=counter_mam_c+1;    end
        if Tmam_e>Th_MAM_e;      counter_mam_e=counter_mam_e+1;    end
        if Tmam_l>Th_MAM_l;      counter_mam_l=counter_mam_l+1;    end
        if Tmam_p>Th_MAM_p;      counter_mam_p=counter_mam_p+1;    end
        if Tmam_ro>Th_MAM_ro;      counter_mam_ro=counter_mam_ro+1;    end
        if Tglrtmam>Th_GLRTMAM;      counter_glrtmam=counter_glrtmam+1;    end
    end
    Pd_SCM_mc(m)=counter_scm/count_Pd;           counter_scm=0;
    Pd_NSCM_mc(m)=counter_nscm/count_Pd;        counter_nscm=0;
    Pd_MAM_r_mc(m)=counter_mam_r/count_Pd;       counter_mam_r=0; 
%     Pd_MAM_c_mc(m)=counter_mam_c/count_Pd;       counter_mam_c=0;
    Pd_MAM_e_mc(m)=counter_mam_e/count_Pd;       counter_mam_e=0;
    Pd_MAM_l_mc(m)=counter_mam_l/count_Pd;       counter_mam_l=0;
    Pd_MAM_p_mc(m)=counter_mam_p/count_Pd;       counter_mam_p=0;
    Pd_MAM_ro_mc(m)=counter_mam_ro/count_Pd;       counter_mam_ro=0;
    Pd_GLRTMAM_mc(m)=counter_glrtmam/count_Pd;       counter_glrtmam=0; 
end
close(h)
toc
figure(1);
hold on
plot(SNRout,Pd_SCM_mc,'b-+','linewidth',2)
plot(SNRout,Pd_NSCM_mc,'k.-','linewidth',2)
plot(SNRout,Pd_MAM_r_mc,'g-s','linewidth',2)
% plot(SNRout,Pd_MAM_c_mc,'k-s','linewidth',2)
plot(SNRout,Pd_MAM_e_mc,'m-s','linewidth',2)
plot(SNRout,Pd_MAM_l_mc,'r-s','linewidth',2)
plot(SNRout,Pd_MAM_p_mc,'c-s','linewidth',2)
plot(SNRout,Pd_MAM_ro_mc,'b-s','linewidth',2)
plot(SNRout,Pd_GLRTMAM_mc,'y-o','linewidth',2)
h_leg = legend('AMF','ANMF','ReimanDistance',...
    'Euclidean Distance','Log Euclidean Distance',...
    'Power-Euclidean distance','Root Euclidean Distance','GLRTMAM');
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
grid on
% str=['Pd_MAM_IPIX_',num2str(n),'K','_',cdfFile_t,'.mat'];
% save(str);