%%%ʵ��һ������ɫ���ص�GLRT�����
%%%%�û����Ĺ�ʽ
%%8ά�ĵ���ʸ���������õ���GL3
%%��˹�����:GLC��SCMN��CC��SCMN��KGLRT��SCM��NSCM��NSCM
%%%�Ǹ�˹�����:��Ҫ�ó�N�ġ�
%%%cf��ÿһ��ʵ�飬R_KAҪ�̶����䡣
clc
clear 
close all
%%%%��������
n = 1; %����������
str_train = 'p';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
sigma_t = 1;
rou = 0.95;  %%Э����������ɵĳ�������
%%%Pd_CLGLRT_2Kmu1lambda3s0.1o1_p��2K��ѵ����Ԫ��Ŀ��mu��lambda��s��ʧ���������
%%o1:opt=1��p��IG�����ϸ�˹
%%%%�����������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
SNRout=-5:1:25; % ���SNR
cos2=0.9;
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
rouR = zeros(N,N);  %%��ʵ���Ӳ�Э����
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% ϵͳ����ʸ��
for i=1:N
    for j=1:N
        rouR(i,j)=rou^abs(i-j);%*exp(1j*2*pi*abs(i-j)*theta_sig);
    end
end
rouR_abs=abs(rouR);
rouR_half=rouR^0.5;
irouR=inv(rouR);
% t = normrnd(1,sigma_t,N,1);%%0~0.5%%ʧ������
% R_KA = rouR.*(t*t');
load R_KA_g_0.1.mat
iR_KA = inv(R_KA);
%%%%����ʸ������
[UU,SS,VV]=svd(irouR*s);
s_v=UU(:,2); %%%%%% ��vt�ڰ׻��ռ�����������s^H*iR*s_v==0
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
%%%%%��ʽ��ʼ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%���޼���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = waitbar(0,'Please wait...');
tic
parfor i = 1:MonteCarloPfa
    warning('off')
%     warning('query')
%     waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
%%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
    Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % �����źŽ������Ӳ�������
    %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%
    R_SCM = (fun_SCM(Train));
    iR_SCM = inv(R_SCM);
    
    R_SCMN = (fun_SCMN(Train));
    iR_SCMN = inv(R_SCMN);
    
    R_NSCM = fun_NSCM(Train);
    iR_NSCM = inv(R_NSCM);
    
    R_NSCMN = fun_NSCMN(Train);
    iR_NSCMN = inv(R_NSCMN);
    
    R_CC = fun_CC(Train,R_SCMN,R_KA);
    iR_CC = inv(R_CC);
    %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% KGLRT
    Tglrt(i) = fun_1SGLRT(R_SCM,x0,s,mu);
    %%%%%% KGLRTCC
    Tglrtcc(i) = fun_1SGLRT(R_CC,x0,s,mu);
    %%%%%% KGLRTNSCM
    Tglrtnscm(i) = fun_1SGLRT(R_NSCM,x0,s,mu);
    %%%%%% CLGLRT
%     Tclglrt(i) = fun_CLGLRT2(lambda,mu,R_KA,R_SCMN,x0,s);
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
%%%%%%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter_glrt=0;
counter_clglrt=0;
counter_glrtcc=0;
counter_glrtnscm=0;

Pd_KGLRT_mc = zeros(1,length(SNRout));
Pd_CLGLRT_mc = zeros(1,length(SNRout));
Pd_KGLRTCC_mc = zeros(1,length(SNRout));
Pd_KGLRTNSCM_mc = zeros(1,length(SNRout));
% alpha=sqrt(SNRnum/abs(s_real'*irouR*s_real)); % ����SNR=|alpha|^2*s'*R^(-1)*s���|alpha|
alpha=sqrt(SNRnum/abs(s'*irouR*s)); % ����SNR=|alpha|^2*s'*R^(-1)*s���|alpha|
h = waitbar(0,'Please wait...');
tic
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
    parfor i=1:MonteCarloPd
        warning('off')
%       waitbar(((m-1)*MonteCarloPd+i)/length(SNRout)/MonteCarloPd,h,sprintf([num2str(((m-1)*MonteCarloPd+i)/length(SNRout)/MonteCarloPd*100),'%%']));
        %%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
        x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % �����źŽ������Ӳ�������
        %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%
        R_SCM = (fun_SCM(Train));
        iR_SCM = inv(R_SCM);
        
        R_SCMN = (fun_SCMN(Train));
        iR_SCMN = inv(R_SCMN);
        
        R_NSCM = (fun_NSCM(Train));
        iR_NSCM = inv(R_NSCM);
        
        R_NSCMN = fun_NSCMN(Train);
        iR_NSCMN = inv(R_NSCMN);
        
        R_CC = fun_CC(Train,R_SCMN,R_KA);
        iR_CC = inv(R_CC);
        x0=alpha(m)*s+x0;%+pp;    %%%%%%%  ��Ҫ  %%%%%%%%%%%%%
        %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%% 1SGLRT
        Tglrt = fun_1SGLRT(R_SCM,x0,s,mu);
        %%%%%% 1SGLRTCC
        Tglrtcc = fun_1SGLRT(R_CC,x0,s,mu);
        %%%%%% 1SGLRTNSCM
        Tglrtnscm = fun_1SGLRT(R_NSCM,x0,s,mu);
        %%%%%% GLCGLRT
%         Tclglrt = fun_CLGLRT2(lambda,mu,R_KA,R_SCMN,x0,s);
        Tclglrt = fun_CLGLRT3(lambda,mu,R_KA,R_SCMN,x0,s);
        %%%�ж�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
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
figure(1);
hold on
plot(SNRout,Pd_CLGLRT_mc,'k.-','linewidth',2)
plot(SNRout,Pd_KGLRTCC_mc,'g-s','linewidth',2)
plot(SNRout,Pd_KGLRT_mc,'b-+','linewidth',2)
plot(SNRout,Pd_KGLRTNSCM_mc,'c-o','linewidth',2)
h_leg = legend('CLGLRT','KGLRTCC','KGLRT','KGLRTNSCM');
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
grid on
% % str=['Pd_CLGLRT_',num2str(N),'N','_',num2str(n),'K','mu',num2str(mu),'lambda',num2str(lambda),'s',num2str(sigma_t),'o',num2str(opt_train),'_',str_train,'.mat'];
str=['Pd_CLGLRT4_3_',num2str(n),'K','mu',num2str(mu),...
     'lambda', num2str(lambda),'s',num2str(sigma_t),...
     'o',num2str(opt_train),'_',str_train,'.mat'];
save(str,'lambda','mu','sigma_t','SNRout','Pd_CLGLRT_mc','Pd_KGLRT_mc',...
         'Pd_KGLRTCC_mc','Pd_KGLRTNSCM_mc');
