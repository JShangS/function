%%%实现一个基于色加载的GLRT检测器
clc
clear 
close all
%%%%假设参数设置
Na = 4;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
SNRout=0:1:40; % 输出SNR
cos2=0.9;
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
rou = 0.95;  %%协方差矩阵生成的迟滞因子
rouR = zeros(N,N);  %%真实的杂波协方差
n = 2;
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
for i=1:N
    for j=1:N
        rouR(i,j)=rou^abs(i-j);%*exp(1j*2*pi*abs(i-j)*theta_sig)
    end
end
irouR=inv(rouR);
rouR_abs=abs(rouR);
sigma_t = 0.1;
t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
R_KA = rouR.*(t*t');
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
% figure; plot(weight,cos2_tmpt);
%%%%%正式开始%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamda = 1.2;
mu = 1;
%%%门限计算
 h = waitbar(0,'Please wait...');
for i = 1:MonteCarloPfa
    waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
%     Train = fun_TrainData_gauss(N,L,rouR);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
%     x0 = fun_TrainData_gauss(N,1,rouR); % 接收信号仅包括杂波和噪声
    Train = fun_TrainData_IGCC(N,L,rouR,lamda,mu,2);%%产生的训练数据,协方差矩阵为rouR的逆gamma纹理复合高斯杂波
    x0 = fun_TrainData_IGCC(N,1,rouR,lamda,mu,2); % 接收信号仅包括杂波和噪声
%     Train = fun_TrainData_K(N,L,rouR,mu);%%产生的训练数据,协方差矩阵为rouR的复合高斯杂波——K分布
%     x0 = fun_TrainData_K(N,1,rouR,mu); % 接收信号仅包括杂波和噪声
    %%%%协方差估计
    R_SCM = (fun_SCM(Train));
    iR_SCM = inv(R_SCM);
%     R_NSCM = fun_NSCM(Train);
%     iR_NSCM = inv(R_NSCM);
    R_CC = fun_CC(Train,R_SCM,R_KA);
    iR_CC = inv(R_CC);
%     R_ICL1 = (fun_ICL(s,x0,inv(R_KA),inv(R_SCM),lamda,mu,1));
%     R_ICL0 = (fun_ICL(s,x0,inv(R_KA),inv(R_SCM),lamda,mu,0));
%     iR_ICL1 = inv(R_ICL1);
%     iR_ICL0 = inv(R_ICL0);
    %%%检测器
    Tamf(i) = abs(s'*iR_SCM*x0)^2/abs(s'*iR_SCM*s);     %%%%%% AMF或者wald
    tmp=abs(x0'*iR_SCM*x0);
%     Tamf(i) = abs(s'*iR_NSCM*x0)^2/abs(s'*iR_NSCM*s);     %%%%%% AMF或者wald
%     tmp=abs(x0'*iR_NSCM*x0);
    Tamf_CC(i) = abs(s'*iR_CC*x0)^2/abs(s'*iR_CC*s);     %%%%%% AMF_CC或者wald
    tmp_CC=abs(x0'*iR_CC*x0);
    Tglrt(i) = Tamf(i)/(1+tmp);                   %%%%%% KGLRT
    Tglrt_CC(i) = Tamf_CC(i)/(1+tmp_CC);         %%%%%% KGLRT_CC
    Tace(i)=Tamf(i)/tmp;                        %%%%%% ACE
    Tabort(i)=(1+Tamf(i))/(2+tmp);              %%%%%% ABORT  % eq.(16) 检测统计量
    Twabort(i)=1/(1+tmp)/(1-Tglrt(i))^2;        %%%%%% ABORT  % 见会议论文中的eq.(18)
    Tace_bar=Tace(i)/(1-Tace(i));
    Tprao(i)=Tglrt(i)^2/(Tamf(i)*(1-Tglrt(i))); %%%%%% DMRao
    Tdnamf(i)=Tace_bar/tmp;                     %%%%%% DNAMF  % eq.(24) 检测统计量
    Taed(i)=tmp;                                %%%%%% 能量检测器 
    %%%%%% CLGLRT
    Tclglrt(i) = fun_CLGRT(lamda,mu,R_KA,R_SCM,x0,s);
%     a = (s'*iR_ICL1*x0)/(s'*iR_ICL1*s);
%     tmp1 = det(iR_ICL1)*((x0 - a*s)'*iR_ICL1*(x0 - a*s)+1/mu)^(-lamda-N);
%     tmp2 = det(iR_ICL0)*(x0'*iR_ICL0*x0+1/mu)^(-lamda-N);
%     Tclglrt(i) =  abs(tmp1/tmp2);%%%%%% 色加载的GLRT
end
close(h)
TAMF=sort(Tamf,'descend');
TACE=sort(Tace,'descend');
TKGLRT=sort(Tglrt,'descend');
TABORT=sort(Tabort,'descend');
TWABORT=sort(Twabort,'descend');
TDMRao=sort(Tprao,'descend');
TDNAMF=sort(Tdnamf,'descend');
TAED=sort(Taed,'descend');
TCLGLRT=sort(Tclglrt,'descend');
TKGLRTCC = sort(Tglrt_CC,'descend');
TAMFCC=sort(Tamf_CC,'descend');

Th_AMF=(TAMF(floor(MonteCarloPfa*PFA-1))+TAMF(floor(MonteCarloPfa*PFA)))/2;
Th_ACE=(TACE(floor(MonteCarloPfa*PFA-1))+TACE(floor(MonteCarloPfa*PFA)))/2;
Th_KGLRT=(TKGLRT(floor(MonteCarloPfa*PFA-1))+TKGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_ABORT=(TABORT(floor(MonteCarloPfa*PFA-1))+TABORT(floor(MonteCarloPfa*PFA)))/2;
Th_WABORT=(TWABORT(floor(MonteCarloPfa*PFA-1))+TWABORT(floor(MonteCarloPfa*PFA)))/2;
Th_DMRao=(TDMRao(floor(MonteCarloPfa*PFA-1))+TDMRao(floor(MonteCarloPfa*PFA)))/2;
Th_DNAMF=(TDNAMF(floor(MonteCarloPfa*PFA-1))+TDNAMF(floor(MonteCarloPfa*PFA)))/2;
Th_AED=(TAED(floor(MonteCarloPfa*PFA-1))+TAED(floor(MonteCarloPfa*PFA)))/2;
Th_CLGLRT=(TCLGLRT(floor(MonteCarloPfa*PFA-1))+TCLGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_KGLRTCC=(TKGLRTCC(floor(MonteCarloPfa*PFA-1))+TKGLRTCC(floor(MonteCarloPfa*PFA)))/2;
Th_AMFCC=(TAMFCC(floor(MonteCarloPfa*PFA-1))+TAMFCC(floor(MonteCarloPfa*PFA)))/2;
%%%%检测曲线%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter_amf=0;
counter_ace=0;
counter_glrt=0;
counter_abort=0;
counter_wabort=0;
counter_prao=0;
counter_dnamf=0;
counter_aed=0;
counter_clglrt=0;
counter_amfcc=0;
counter_glrtcc=0;
Pd_AMF_mc = zeros(1,length(SNRout));
Pd_KGLRT_mc = zeros(1,length(SNRout));
Pd_ACE_mc = zeros(1,length(SNRout));
Pd_ABORT_mc = zeros(1,length(SNRout));
Pd_WABORT_mc = zeros(1,length(SNRout));
Pd_DMRao_mc = zeros(1,length(SNRout));
Pd_DNAMF_mc = zeros(1,length(SNRout));
Pd_AED_mc = zeros(1,length(SNRout));
Pd_CLGLRT_mc = zeros(1,length(SNRout));
Pd_KGLRTCC_mc = zeros(1,length(SNRout));
Pd_AMFCC_mc = zeros(1,length(SNRout));
alpha=sqrt(SNRnum/abs(s_real'*irouR*s_real)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
h = waitbar(0,'Please wait...');
for m=1:length(SNRout)
    m
    for i=1:MonteCarloPd 
        waitbar(((m-1)*MonteCarloPd+i)/length(SNRout)/MonteCarloPd,h,sprintf([num2str(((m-1)*MonteCarloPd+i)/length(SNRout)/MonteCarloPd*100),'%%']));
%         Train = fun_TrainData_gauss(N,L,rouR);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
%         x0 = fun_TrainData_gauss(N,1,rouR); % 接收信号仅包括杂波和噪声
        Train = fun_TrainData_IGCC(N,L,rouR,lamda,mu,2);%%产生的训练数据,协方差矩阵为rouR的逆gamma纹理复合高斯杂波
        x0 = fun_TrainData_IGCC(N,1,rouR,lamda,mu,2); % 接收信号仅包括杂波和噪声
%         Train = fun_TrainData_K(N,L,rouR,mu);%%产生的训练数据,协方差矩阵为rouR的复合高斯杂波——K分布
%         x0 = fun_TrainData_K(N,1,rouR,mu); % 接收信号仅包括杂波和噪声
        %%%%协方差估计
        R_SCM = (fun_SCM(Train));
        iR_SCM = inv(R_SCM);
%         R_NSCM = (fun_NSCM(Train));
%         iR_NSCM = inv(R_NSCM);
        R_CC = fun_CC(Train,R_SCM,R_KA);
        iR_CC = inv(R_CC);
        x0=alpha(m)*s_real+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
%         R_ICL1 = (fun_ICL(s,x0,inv(R_KA),inv(R_SCM),lamda,mu,1));
%         R_ICL0 = (fun_ICL(s,x0,inv(R_KA),inv(R_SCM),lamda,mu,0));
%         iR_ICL1 = inv(R_ICL1);
%         iR_ICL0 = inv(R_ICL0);
        %%%检测器
        Tamf = abs(s'*iR_SCM*x0)^2/abs(s'*iR_SCM*s);     %%%%%% AMF或者wald
        tmp=abs(x0'*iR_SCM*x0);
%         Tamf = abs(s'*iR_NSCM*x0)^2/abs(s'*iR_NSCM*s);     %%%%%% AMF或者wald
%         tmp=abs(x0'*iR_NSCM*x0);
        Tamf_CC = abs(s'*iR_CC*x0)^2/abs(s'*iR_CC*s);     %%%%%% AMF或者wald
        tmp_CC=abs(x0'*iR_CC*x0);
        Tglrt = Tamf/(1+tmp);                   %%%%%% KGLRT
        Tglrt_CC = Tamf_CC/(1+tmp_CC);   
        Tace=Tamf/tmp;                        %%%%%% ACE
        Tabort=(1+Tamf)/(2+tmp);              %%%%%% ABORT  % eq.(16) 检测统计量
        Twabort=1/(1+tmp)/(1-Tglrt)^2;        %%%%%% ABORT  % 见会议论文中的eq.(18)
        Tace_bar=Tace/(1-Tace);
        Tprao=Tglrt^2/(Tamf*(1-Tglrt));       %%%%%% DMRao
        Tdnamf=Tace_bar/tmp;                  %%%%%% DNAMF  % eq.(24) 检测统计量
        Taed=tmp;                             %%%%%% 能量检测器  
        %%%%%% CLGLRT
          Tclglrt = fun_CLGRT(lamda,mu,R_KA,R_SCM,x0,s);
%         a = (s'*iR_ICL1*x0)/(s'*iR_ICL1*s);
%         tmp1 = det(iR_ICL1)*((x0 - a*s)'*iR_ICL1*(x0 - a*s)+1/mu)^(-lamda-N);
%         tmp2 = det(iR_ICL0)*(x0'*iR_ICL0*x0+1/mu)^(-lamda-N);
%         Tclglrt =  abs(tmp1/tmp2);%%%%%% 色加载的GLRT
        if Tamf>Th_AMF;         counter_amf=counter_amf+1;          end            
        if Tglrt>Th_KGLRT;      counter_glrt=counter_glrt+1;        end                
        if Tace>Th_ACE;         counter_ace=counter_ace+1;          end          
        if Tabort>Th_ABORT;     counter_abort=counter_abort+1;      end            
        if Twabort>Th_WABORT;   counter_wabort=counter_wabort+1;    end        
        if Tprao>Th_DMRao;      counter_prao=counter_prao+1;        end          
        if Tdnamf>Th_DNAMF;     counter_dnamf=counter_dnamf+1;      end          
        if Taed>Th_AED;         counter_aed=counter_aed+1;          end            
        if Tclglrt>Th_CLGLRT;      counter_clglrt=counter_clglrt+1;        end
        if Tglrt_CC>Th_KGLRTCC;      counter_glrtcc=counter_glrtcc+1;        end 
        if Tamf_CC>Th_AMFCC;      counter_amfcc=counter_amfcc+1;        end 
    end
    Pd_AMF_mc(m)=counter_amf/MonteCarloPd;          counter_amf=0;
    Pd_KGLRT_mc(m)=counter_glrt/MonteCarloPd;        counter_glrt=0;
    Pd_ACE_mc(m)=counter_ace/MonteCarloPd;          counter_ace=0;
    Pd_ABORT_mc(m)=counter_abort/MonteCarloPd;      counter_abort=0;
    Pd_WABORT_mc(m)=counter_wabort/MonteCarloPd;    counter_wabort=0;
    Pd_DMRao_mc(m)=counter_prao/MonteCarloPd;        counter_prao=0;
    Pd_DNAMF_mc(m)=counter_dnamf/MonteCarloPd;      counter_dnamf=0;
    Pd_AED_mc(m)=counter_aed/MonteCarloPd;          counter_aed=0; 
    Pd_CLGLRT_mc(m)=counter_clglrt/MonteCarloPd;        counter_clglrt=0;
    Pd_KGLRTCC_mc(m)=counter_glrtcc/MonteCarloPd;        counter_glrtcc=0;
    Pd_AMFCC_mc(m)=counter_amfcc/MonteCarloPd;        counter_amfcc=0;
end
close(h)
figure(2);
hold on
plot(SNRout,Pd_KGLRT_mc,'b-+','linewidth',2)
plot(SNRout,Pd_AMF_mc,'g.-')
plot(SNRout,Pd_ACE_mc,'r-x','linewidth',2)
plot(SNRout,Pd_ABORT_mc,'c-*','linewidth',2)
plot(SNRout,Pd_WABORT_mc,'m-P','linewidth',2)
plot(SNRout,Pd_DMRao_mc,'r-o','linewidth',2)
plot(SNRout,Pd_AED_mc,'r-d','linewidth',2)
plot(SNRout,Pd_DNAMF_mc,'g-s','linewidth',2);
plot(SNRout,Pd_CLGLRT_mc,'k.-')
plot(SNRout,Pd_KGLRTCC_mc,'b-s')
plot(SNRout,Pd_AMFCC_mc,'g-o')
legend('KGLRT','AMF/DMwald','ACE','ABORT','WABORT','DMRao','AED','DNAMF','CLGLRT','KGLRTCC','AMFCC')
% legend({'KGLRT','AMF/DMwald','DMRao'},'FontSize',20)
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
grid on
str=['Pd_CLGLRT_',num2str(n),'K','_mu',num2str(mu),'_lamda',num2str(lamda),'_sigma',num2str(sigma_t),'.mat'];
save(str,'lamda','mu','SNRout','Pd_ABORT_mc','Pd_ACE_mc','Pd_AED_mc','Pd_AMF_mc',...
    'Pd_CLGLRT_mc','Pd_DMRao_mc','Pd_DNAMF_mc','Pd_KGLRT_mc','Pd_WABORT_mc',...
    'Pd_KGLRTCC_mc', 'Pd_AMFCC_mc');
