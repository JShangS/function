%%%GLC-GLRT的ROC曲线
clc
clear 
close all
% %%%%参数设置
n = 2; %几倍的样本
str_train = 'p';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
sigma_t = 0.1;
%%%Pd_CLGLRT_2Kmu1lambda3s0.1o1_p：2K：训练单元数目，mu，lambda，s：失配向量方差，
%%o1:opt=1，p：IG纹理复合高斯
%%%%假设参数设置
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
SNRout=[-5,5,15]; % 输出SNR
cos2=0.9;
PFA=[1e-4,1e-3:1e-2:1e-1+1e-2];% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=round(1./PFA*100);
L_Pfa = length(MonteCarloPfa);
MonteCarloPd=1e4;
rou = 0.95;  %%协方差矩阵生成的迟滞因子
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
R_KA = zeros(size(rouR));
tic
t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
R_KA = R_KA+rouR.*(t*t')/10000;
iR_KA = inv(R_KA);
toc
% R_KA_inv = inv(R_KA);
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
% %%%%%正式开始%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%门限计算%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0,'Please wait...');
for i_Pfa = 1:L_Pfa
    waitbar((i_Pfa/L_Pfa),h,sprintf([num2str(i_Pfa/L_Pfa*100),'%%']));
    Tamf = zeros(MonteCarloPfa(i_Pfa),1);
    Tamfcc = zeros(MonteCarloPfa(i_Pfa),1);
    Tamfnscm = zeros(MonteCarloPfa(i_Pfa),1);
    Tglrt = zeros(MonteCarloPfa(i_Pfa),1);
    Tglrtcc = zeros(MonteCarloPfa(i_Pfa),1);
    Tglrtnscm = zeros(MonteCarloPfa(i_Pfa),1);
    Tclglrt = zeros(MonteCarloPfa(i_Pfa),1);
    parfor i = 1:MonteCarloPfa(i_Pfa)
%     waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
%%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
    Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
    %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%
    R_SCM = (fun_SCM(Train));
    iR_SCM = inv(R_SCM);
    
    R_SCMN = (fun_SCMN(Train));
    iR_SCMN = inv(R_SCMN);
    
    R_NSCM = fun_NSCM(Train);
    iR_NSCM = inv(R_NSCM);
    
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
    Tclglrt(i) = fun_CLGLRT3(lambda,mu,R_KA,R_SCMN,x0,s);
    end
    % close(h)
    TKGLRT=sort(Tglrt,'descend');
    TCLGLRT=sort(Tclglrt,'descend');
    TKGLRTCC=sort(Tglrtcc,'descend');
    TKGLRTNSCM=sort(Tglrtnscm,'descend');

    Th_KGLRT(i_Pfa)=(TKGLRT(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa)-1))+TKGLRT(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa))))/2;
    Th_CLGLRT(i_Pfa)=(TCLGLRT(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa)-1))+TCLGLRT(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa))))/2;
    Th_KGLRTCC(i_Pfa)=(TKGLRTCC(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa)-1))+TKGLRTCC(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa))))/2;
    Th_KGLRTNSCM(i_Pfa)=(TKGLRTNSCM(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa)-1))+TKGLRTNSCM(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa))))/2;
end
close(h)
save('ThPfa2','PFA','Th_KGLRT','Th_CLGLRT','Th_KGLRTCC','Th_KGLRTNSCM','L',...
    'str_train','lambda','mu');
% %%%%%%%%%%%%%%%%%%%%%检测概率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load ThPfa.mat
Pd_KGLRT_Mlti_mc = zeros(L_Pfa,length(SNRout));
Pd_CLGLRT_Mlti_mc = zeros(L_Pfa,length(SNRout));
Pd_KGLRTCC_Mlti_mc = zeros(L_Pfa,length(SNRout));
Pd_KGLRTNSCM_Mlti_mc = zeros(L_Pfa,length(SNRout));

counter_glrt=0;
counter_clglrt=0;
counter_glrtcc=0;
counter_amfcc=0;
counter_glrtnscm=0;
alpha=sqrt(SNRnum/abs(s_real'*irouR*s_real)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
h = waitbar(0,'Please wait...');
tic
L_SNRout = length(SNRout);
for i_Pfa = 1:L_Pfa %%虚警
   for m=1:L_SNRout %%信噪比
       waitbar(((i_Pfa-1)*L_SNRout+m)/(L_SNRout*L_Pfa),h,sprintf([num2str(((i_Pfa-1)*L_SNRout+m)/(L_SNRout*L_Pfa)*100),'%%']));
        parfor i=1:MonteCarloPd %%%MC检测概率
    %         waitbar(((m-1)*MonteCarloPd+i)/length(SNRout)/MonteCarloPd,h,sprintf([num2str(((m-1)*MonteCarloPd+i)/length(SNRout)/MonteCarloPd*100),'%%']));
            %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
            Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
            x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
            %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%
            R_SCM = (fun_SCM(Train));
            iR_SCM = inv(R_SCM);
            R_SCMN = (fun_SCMN(Train));
            iR_SCMN = inv(R_SCMN);
            R_NSCM = (fun_NSCM(Train));
            iR_NSCM = inv(R_NSCM);
            R_CC = fun_CC(Train,R_SCMN,R_KA);
            iR_CC = inv(R_CC);
            x0=alpha(m)*s_real+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
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
            Tclglrt = fun_CLGLRT3(lambda,mu,R_KA,R_SCMN,x0,s);
            %%%判断%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            if Tglrt>Th_KGLRT(i_Pfa);          counter_glrt=counter_glrt+1;        end                  
            if Tclglrt>Th_CLGLRT(i_Pfa);       counter_clglrt=counter_clglrt+1;    end   
            if Tglrtcc>Th_KGLRTCC(i_Pfa);      counter_glrtcc=counter_glrtcc+1;    end
            if Tglrtnscm>Th_KGLRTNSCM(i_Pfa);  counter_glrtnscm=counter_glrtnscm+1;      end
        end
        Pd_KGLRT_Mlti_mc(i_Pfa,m)=counter_glrt/MonteCarloPd;           counter_glrt=0;
        Pd_CLGLRT_Mlti_mc(i_Pfa,m)=counter_clglrt/MonteCarloPd;        counter_clglrt=0;
        Pd_KGLRTCC_Mlti_mc(i_Pfa,m)=counter_glrtcc/MonteCarloPd;       counter_glrtcc=0;
        Pd_KGLRTNSCM_Mlti_mc(i_Pfa,m)=counter_glrtnscm/MonteCarloPd;    counter_glrtnscm=0;
    end
end

close(h)
toc
figure(2);
hold on
% plot(PFA,Pd_CLGLRT_Mlti_mc(:,1),'ks','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_KGLRTCC_Mlti_mc(:,1),'ko','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_KGLRT_Mlti_mc(:,1),'k>','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_KGLRTNSCM_Mlti_mc(:,1),'k*','linewidth',2,'MarkerSize',15)
%%5dB
plot(PFA,Pd_CLGLRT_Mlti_mc(:,1),'k-s','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRTCC_Mlti_mc(:,1),'k-o','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRT_Mlti_mc(:,1),'k->','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRTNSCM_Mlti_mc(:,1),'k-*','linewidth',2,'MarkerSize',15)
%%10dB
plot(PFA,Pd_CLGLRT_Mlti_mc(:,2),'r-s','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRTCC_Mlti_mc(:,2),'r-o','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRT_Mlti_mc(:,2),'r->','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRTNSCM_Mlti_mc(:,2),'r-*','linewidth',2,'MarkerSize',15)
%%15dB
plot(PFA,Pd_CLGLRT_Mlti_mc(:,3),'g-s','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRTCC_Mlti_mc(:,3),'g-o','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRT_Mlti_mc(:,3),'g->','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRTNSCM_Mlti_mc(:,3),'g-*','linewidth',2,'MarkerSize',15)
h_leg = legend('GLC-GLRT','GLRT with CC','GLRT with SCM','GLRT with NSCM');
xlabel('Pfa','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
grid on
box on
str=['Pd_CLGLRT2_ROC3',num2str(n),'K','mu',num2str(mu),...
     'lambda', num2str(lambda),'s',num2str(sigma_t),...
     'o',num2str(opt_train),'_',str_train,'.mat'];
save(str,'lambda','mu','sigma_t','SNRout','PFA','Pd_CLGLRT_Mlti_mc','Pd_KGLRT_Mlti_mc',...
    'Pd_KGLRTCC_Mlti_mc','Pd_KGLRTNSCM_Mlti_mc');
