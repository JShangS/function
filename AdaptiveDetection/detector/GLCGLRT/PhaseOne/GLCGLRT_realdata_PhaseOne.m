%%%实现一个基于色加载的GLRT检测器
%%%%MIT Lincoln Laboratory Phase One radar 1985年数据
%%%%《CFAR Behavior of Adaptive Detectors: An Experimental Analysis》
clc
clear 
close all
load Real_data_PhaseOneRadar_range30.mat
%%%%参数设置
n = 2; %几倍的样本
sigma_t = 0;
lambda = 3;
mu = 1;
str_train = 'g';
opt_train = 1;
%%%%参数设置
N = 8;
SNRout=-5:1:15; % 输出SNR
cos2=0.9;
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
% rou = 0.95;  %%协方差矩阵生成的迟滞因子
rouR = R_KA;  %%真实的杂波协方差
irouR = inv(rouR);
L=round(n*N); 
theta_sig = 0.3;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
R_KA = rouR.*(t*t');
iR_KA = inv(R_KA);
rouR_half=rouR^0.5;
%%%%导向矢量设置
[UU,SS,VV]=svd(irouR*s);
s_v=UU(:,2); %%%%%% 与vt在白化空间正交，即：s^H*iR*s_v==0
weight=linspace(0,1,300);
for i=1:length(weight)
    s_tmpt=weight(i)*s+(1-weight(i))*s_v;
    cos2_tmpt(i)=abs(s_tmpt'*irouR*s).^2/abs(s_tmpt'*irouR*s_tmpt*s'*irouR*s);
end
[Min, Index]=min(abs(cos2-cos2_tmpt));
Weight=weight(Index);
s_real=Weight*s+(1-Weight)*s_v;
% figure;plot(abs(s_real))
% figure; plot(weight,cos2_tmpt);
%%%%%正式开始%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%门限计算%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = waitbar(0,'Please wait...');
tic
M = 30720;
parfor i = 1:MonteCarloPfa
%     waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
%%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%   
%     Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
%     x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
    index_t1 = ceil(rand()*(M-10));
    Train1 = Zhh(index_t1:index_t1+7,Range-L/2:Range-1);
    Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2);
    Train = [Train1,Train2];%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = Zhh(index_t1:index_t1+7,Range) ; % 接收信号仅包括杂波和噪声
    %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%
    R_SCM = (fun_SCM(Train));
    iR_SCM = inv(R_SCM);
    
    R_SCMN = (fun_SCMN(Train));
    iR_SCMN = inv(R_SCMN);
    
    R_NSCM = fun_NSCM(Train);
    iR_NSCM = inv(R_NSCM);
    
    R_NSCMN = fun_NSCM(Train);
    iR_NSCMN = inv(R_NSCMN);
    
    R_CC = fun_CC(Train,R_SCMN,R_KA);
    iR_CC = inv(R_CC);
    %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% AMF或者wald
    Tamf(i) = abs(s'*iR_SCM*x0)^2/abs(s'*iR_SCM*s);    
    tmp=abs(x0'*iR_SCM*x0);
    %%%%%% AMFCC或者wald_CC
    Tamfcc(i) = abs(s'*iR_CC*x0)^2/abs(s'*iR_CC*s);     
    tmpcc=abs(x0'*iR_CC*x0);
    %%%%%%%%%%% AMFNSCM
    Tamfnscm(i) = abs(s'*iR_NSCM*x0)^2/abs(s'*iR_NSCM*s);     
    tmpnscm=abs(x0'*iR_NSCM*x0);
    %%%%%% KGLRT
    Tglrt(i) = Tamf(i)/(1+tmp);     
    %%%%%% KGLRTCC
    Tglrtcc(i) = Tamfcc(i)/(1+tmpcc);
    %%%%%% KGLRTNSCM
    Tglrtnscm(i) = Tamfnscm(i)/(1+tmpnscm);
    %%%%%% CLGLRT
%     Tclglrt(i) = fun_CLGLRT2(R_KA,R_SCM,x0,s);
    Tclglrt(i) = fun_CLGLRT3(lambda,mu,R_KA,R_SCMN,x0,s);
end
toc
% close(h)
TKGLRT=sort(Tglrt,'descend');
TCLGLRT=sort(Tclglrt,'descend');
TKGLRTCC=sort(Tglrtcc,'descend');
TKGLRTNSCM=sort(Tglrtnscm,'descend');

Th_KGLRT=(TKGLRT(floor(MonteCarloPfa*PFA-1))+TKGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_CLGLRT=(TCLGLRT(floor(MonteCarloPfa*PFA-1))+TCLGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_KGLRTCC=(TKGLRTCC(floor(MonteCarloPfa*PFA-1))+TKGLRTCC(floor(MonteCarloPfa*PFA)))/2;
Th_KGLRTNSCM=(TKGLRTNSCM(floor(MonteCarloPfa*PFA-1))+TKGLRTNSCM(floor(MonteCarloPfa*PFA)))/2;
%%%%%%%%%%%%%%%%%%%%%检测概率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

counter_glrt=0;
counter_clglrt=0;
counter_glrtcc=0;
counter_glrtnscm=0;

Pd_KGLRT_mc = zeros(1,length(SNRout));
Pd_CLGLRT_mc = zeros(1,length(SNRout));
Pd_KGLRTCC_mc = zeros(1,length(SNRout));
Pd_KGLRTNSCM_mc = zeros(1,length(SNRout));
alpha=sqrt(SNRnum/abs(s_real'*irouR*s_real)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
% alpha=sqrt(SNRnum/abs(s'*irouR*s));
h = waitbar(0,'Please wait...');
tic
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
    parfor i=1:MonteCarloPd 
%         waitbar(((m-1)*MonteCarloPd+i)/length(SNRout)/MonteCarloPd,h,sprintf([num2str(((m-1)*MonteCarloPd+i)/length(SNRout)/MonteCarloPd*100),'%%']));
        %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
%         Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
%         x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
        index_t1 = ceil(rand()*(M-10));
        Train1 = Zhh(index_t1:index_t1+7,Range-L/2:Range-1);
        Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2);
        Train = [Train1,Train2];%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        x0 = Zhh(index_t1:index_t1+7,Range) ; % 接收信号仅包括杂波和噪声
        x0=alpha(m)*s+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
        %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%
        R_SCM = (fun_SCM(Train));
        iR_SCM = inv(R_SCM);
        
        R_SCMN = (fun_SCMN(Train));
        iR_SCMN = inv(R_SCMN);
        
        R_NSCM = (fun_NSCM(Train));
        iR_NSCM = inv(R_NSCM);
        
        R_NSCMN = fun_NSCM(Train);
        iR_NSCMN = inv(R_NSCMN);
        R_CC = fun_CC(Train,R_SCMN,R_KA);
        iR_CC = inv(R_CC);
        
        %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% AMF或者wald
        Tamf = abs(s'*iR_SCM*x0)^2/abs(s'*iR_SCM*s);   
        tmp=abs(x0'*iR_SCM*x0);
        %%%%%% AMFCC或者wald
        Tamfcc = abs(s'*iR_CC*x0)^2/abs(s'*iR_CC*s);    
        tmpcc = abs(x0'*iR_CC*x0);
        %%%%%% AMF-NSCM或者wald即NAMF
        Tamfnscm = abs(s'*iR_NSCM*x0)^2/abs(s'*iR_NSCM*s);    
        tmpnscm = abs(x0'*iR_NSCM*x0);
        %%%%%% KGLRT
        Tglrt = Tamf/(1+tmp); 
        %%%%%% KGLRTCC
        Tglrtcc = Tamfcc/(1+tmpcc);
        %%%%%% KGLRTNSCM
        Tglrtnscm = Tamfnscm/(1+tmpnscm);       
        %%%%%% CLGLRT
%         Tclglrt = fun_CLGLRT2(R_KA,R_SCM,x0,s);
        Tclglrt = fun_CLGLRT3(lambda,mu,R_KA,R_SCMN,x0,s);
        %%%判断%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        if Tglrt>Th_KGLRT;          counter_glrt=counter_glrt+1;        end                          
        if Tclglrt>Th_CLGLRT;       counter_clglrt=counter_clglrt+1;    end   
        if Tglrtcc>Th_KGLRTCC;      counter_glrtcc=counter_glrtcc+1;    end
        if Tglrtnscm>Th_KGLRTNSCM;  counter_glrtnscm=counter_glrtnscm+1;      end
    end
    Pd_KGLRT_mc(m)=counter_glrt/MonteCarloPd;           counter_glrt=0;
    Pd_CLGLRT_mc(m)=counter_clglrt/MonteCarloPd;        counter_clglrt=0;
    Pd_KGLRTCC_mc(m)=counter_glrtcc/MonteCarloPd;       counter_glrtcc=0;
    Pd_KGLRTNSCM_mc(m)=counter_glrtnscm/MonteCarloPd;    counter_glrtnscm=0;
end
close(h)
toc
figure(2);
hold on

plot(SNRout,Pd_CLGLRT_mc,'k-p','linewidth',2)
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',2)
plot(SNRout,Pd_KGLRT_mc,'b-o','linewidth',2)
plot(SNRout,Pd_KGLRTNSCM_mc,'c-s','linewidth',2)
legend('CLGLRT','KGLRTCC','KGLRT','KGLRTNSCM');
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
grid on

% str=['Pd_CLGLRT4_PhaseOne_',num2str(n),'K','s',num2str(sigma_t),'.mat'];
% save(str,'SNRout','Pd_CLGLRT_mc','Pd_KGLRT_mc','Pd_KGLRTCC_mc','Pd_KGLRTNSCM_mc');
