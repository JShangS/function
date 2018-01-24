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
n = 1; %几倍的样本
sigma_t = 0.0;
% range = 14;
%%参数估计来自《Maximum Likelihood Estimation for
%%%            Compound-Gaussian Clutter with Inverse GammaTexture》

% lambda = 3;
% mu = 1;
lambda = 1.4495;
mu = 1.3820;
str_train = 'g';
opt_train = 1;
%%%%参数设置
% N = 8;
SNRout=-5:1:25; % 输出SNR
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=100/PFA;
MonteCarloPd=1e4;
% rou = 0.95;  %%协方差矩阵生成的迟滞因子
% load rou_19980223_170435.mat
rouR = R_KA;  %%真实的杂波协方差
irouR = inv(rouR);
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
% t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
% load R_KA.mat
% R_KA = rouR.*(t*t');
iR_KA = inv(R_KA);
rouR_half=rouR^0.5;
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
    iR_SCM = inv(R_SCM);
    
    R_SCMN = (fun_SCMN(Train));
    iR_SCMN = inv(R_SCMN);
    
    R_NSCM = fun_NSCM(Train,1);
    iR_NSCM = inv(R_NSCM);
    
    R_NSCMN = fun_NSCMN(Train);
    iR_NSCMN = inv(R_NSCMN);
    
    R_CC = fun_CC(Train,R_SCMN,R_KA);
    iR_CC = inv(R_CC);
    %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% KGLRT
    Tglrt(i) = fun_1SGLRT(R_SCMN,x0,s,mu);
    %%%%%% KGLRTCC
    Tglrtcc(i) = fun_1SGLRT(R_CC,x0,s,mu);
    %%%%%% KGLRTNSCM
    Tglrtnscm(i) = fun_1SGLRT(R_NSCMN,x0,s,mu);
    %%%%%% CLGLRT
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
alpha=sqrt(SNRnum/abs(s'*irouR*s)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
% alpha=sqrt(SNRnum/abs(s'*irouR*s));
h = waitbar(0,'Please wait...');
tic
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
    parfor i=1:MonteCarloPd 
%         waitbar(((m-1)*MonteCarloPd+i)/length(SNRout)/MonteCarloPd,h,sprintf([num2str(((m-1)*MonteCarloPd+i)/length(SNRout)/MonteCarloPd*100),'%%']));
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
        iR_SCM = inv(R_SCM);
        
        R_SCMN = (fun_SCMN(Train));
        iR_SCMN = inv(R_SCMN);
        
        R_NSCM = fun_NSCM(Train,1);
        iR_NSCM = inv(R_NSCM);
        
        R_NSCMN = fun_NSCMN(Train);
        iR_NSCMN = inv(R_NSCMN);
        
        R_CC = fun_CC(Train,R_SCMN,R_KA);
        iR_CC = inv(R_CC);
        %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% KGLRT
        Tglrt = fun_1SGLRT(R_SCMN,x0,s,mu);
        %%%%%% KGLRTCC
        Tglrtcc = fun_1SGLRT(R_CC,x0,s,mu);
        %%%%%% KGLRTNSCM
        Tglrtnscm = fun_1SGLRT(R_NSCMN,x0,s,mu);
        %%%%%% CLGLRT
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
h_leg = legend('GLCLRT','1S-GLRT with CC','1S-GLRT with SCM','1S-GLRT with NSCM');
% legend({'KGLRT','AMF/DMwald','DMRao'},'FontSize',20)
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
grid on
box on
matFile(end-4:end)=[];
% str=['Pd_CLGLRT_',matFile,'_',num2str(n),'K','s',num2str(sigma_t),'_',str_train,'.mat'];
str=['Pd_CLGLRT_2_',matFile,'_',num2str(n),'K','s',num2str(sigma_t),'.mat'];
save(str,'SNRout','Pd_CLGLRT_mc','Pd_KGLRT_mc','Pd_KGLRTCC_mc','Pd_KGLRTNSCM_mc');
