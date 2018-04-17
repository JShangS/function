clc
clear 
close all
%%%%��������
n = 1; %����������
str_train = 'g';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
rou = 0.95;  %%Э����������ɵĳ�������
sigma_t = 0.5;
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
s = exp(-1i*2*pi*nn*theta_sig)'/sqrt(N); %%%%%% ϵͳ����ʸ��
rouR = fun_rho(rou,N,1,0.2);
rouR_abs=abs(rouR);
rouR_half=rouR^0.5;
irouR=inv(rouR);
t = normrnd(1,sigma_t,N,1);%%0~0.5%%ʧ������
R_KA = zeros(size(rouR));
for i = 1:100
    t = normrnd(1,sigma_t,N,1);%%0~0.5%%ʧ������
    R_KA = R_KA + rouR.*(t*t')/100;
end
tic
parfor i = 1:MonteCarloPfa
    warning('off')
%     warning('query')
%     waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
%%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
    Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % �����źŽ������Ӳ�������
    %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%
    R_x0 = abs(fun_SCMN(x0));
     
    R_SCM = (fun_SCM(Train));
    
    R_SCMN = (fun_SCMN(Train));
    
    R_NSCM = (fun_NSCM(Train));
    
    R_NSCMN = (fun_NSCMN(Train));
    
    R_KL = fun_KLPW(R_KA,R_SCMN);
    
    R_CC = fun_CC(Train,R_SCMN,R_KA);
    
    R_H = 0.5 * R_KA + 0.5 * R_SCM;
    %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% ANMF_SCM
    Tanmf_SCM(i) = fun_AMF(R_SCMN,x0,s);
    %%%%%% ANMF_NSCM
    Tanmf_NSCM(i) = fun_ANMF(R_SCMN,x0,s);
    %%%%%% NMF
    Tnmf(i) = fun_ANMF(rouR,x0,s);
    %%%%%% ANMF_KL
    Tanmf_KL(i) = fun_ANMF(R_KL,x0,s);
    %%%%%% ANMF_CC
    Tanmf_CC(i) = fun_ANMF(R_CC,x0,s);
    %%%%%% ANMF_half
    Tanmf_H(i) = fun_ANMF(R_H,x0,s);
end
toc
TANMF_SCM=sort(Tanmf_SCM,'descend');
TANMF_NSCM=sort(Tanmf_NSCM,'descend');
TNMF=sort(Tnmf,'descend');
TANMF_KL=sort(Tanmf_KL,'descend');
TANMF_CC=sort(Tanmf_CC,'descend');
TANMF_H=sort(Tanmf_H,'descend');

Th_SCM=(TANMF_SCM(floor(MonteCarloPfa*PFA-1))+TANMF_SCM(floor(MonteCarloPfa*PFA)))/2;
Th_NSCM=(TANMF_NSCM(floor(MonteCarloPfa*PFA-1))+TANMF_NSCM(floor(MonteCarloPfa*PFA)))/2;
Th_NMF=(TNMF(floor(MonteCarloPfa*PFA-1))+TNMF(floor(MonteCarloPfa*PFA)))/2;

Th_KL=(TANMF_KL(floor(MonteCarloPfa*PFA-1))+TANMF_KL(floor(MonteCarloPfa*PFA)))/2;
Th_CC=(TANMF_CC(floor(MonteCarloPfa*PFA-1))+TANMF_CC(floor(MonteCarloPfa*PFA)))/2;
Th_H=(TANMF_H(floor(MonteCarloPfa*PFA-1))+TANMF_H(floor(MonteCarloPfa*PFA)))/2;

%%%%%%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter_scm=0;
counter_nscm=0;
counter_nmf=0;
counter_kl=0;
counter_cc=0;
counter_h=0;

Pd_SCM_mc = zeros(1,length(SNRout));
Pd_NSCM_mc = zeros(1,length(SNRout));
Pd_NMF_mc = zeros(1,length(SNRout));
Pd_KL_mc = zeros(1,length(SNRout));
Pd_CC_mc = zeros(1,length(SNRout));
Pd_H_mc = zeros(1,length(SNRout));

alpha=sqrt(SNRnum/abs(s'*irouR*s)); % ����SNR=|alpha|^2*s'*R^(-1)*s���|alpha|
h = waitbar(0,'Please wait...');
tic
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
%     MonteCarloPd = 1e4;
    parfor i=1:MonteCarloPd
        %%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
        x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % �����źŽ������Ӳ�������
        %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%
        R_x0 = abs(fun_SCMN(x0));
        
        R_SCM = (fun_SCM(Train));
    
        R_SCMN = (fun_SCMN(Train));
    
        R_NSCM = (fun_NSCM(Train));
    
        R_NSCMN = (fun_NSCMN(Train));
            
        R_KL = fun_KLPW(R_KA,R_SCM);
    
        R_CC = fun_CC(Train,R_SCMN,R_KA);
        
        R_H = 0.5 * R_KA + 0.5 * R_SCM;
        %%%����ź�
        x0=alpha(m)*s+x0;%+pp;    %%%%%%%  ��Ҫ  %%%%%%%%%%%%%
        %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% AMF
        Tscm = fun_AMF(R_SCMN,x0,s);
        %%%%%% ANMF_NSCM
        Tnscm = fun_ANMF(R_SCMN,x0,s);
        %%%%%% NMF
        Tnmf = fun_ANMF(rouR,x0,s);
        %%%%%% ANMF_KL
        Tkl = fun_ANMF(R_KL,x0,s);
        %%%%%% ANMF_CC
        Tcc = fun_ANMF(R_CC,x0,s);
        %%%%%% ANMF_CC
        Thh = fun_ANMF(R_H,x0,s);
        %%%�ж�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        if Tscm>Th_SCM;          counter_scm=counter_scm+1;        end                
        if Tnscm>Th_NSCM;       counter_nscm=counter_nscm+1;    end   
        if Tnmf>Th_NMF;       counter_nmf=counter_nmf+1;    end
        if Tkl>Th_KL;       counter_kl=counter_kl+1;    end
        if Tnmf>Th_CC;       counter_cc=counter_cc+1;    end
        if Thh>Th_H;       counter_h=counter_h+1;    end
    end
    Pd_SCM_mc(m)=counter_scm/MonteCarloPd;           counter_scm=0;
    Pd_NSCM_mc(m)=counter_nscm/MonteCarloPd;        counter_nscm=0;
    Pd_NMF_mc(m)=counter_nmf/MonteCarloPd;        counter_nmf=0;
    Pd_CC_mc(m)=counter_cc/MonteCarloPd;           counter_cc=0;
    Pd_KL_mc(m)=counter_kl/MonteCarloPd;        counter_kl=0;
    Pd_H_mc(m)=counter_h/MonteCarloPd;        counter_h=0;
end
toc
close(h)
figure(1);
hold on
plot(SNRout,Pd_NMF_mc,'k','linewidth',1)
plot(SNRout,Pd_SCM_mc,'k-+','linewidth',1)
plot(SNRout,Pd_NSCM_mc,'k.-','linewidth',1)
plot(SNRout,Pd_KL_mc,'k-*','linewidth',1)
plot(SNRout,Pd_CC_mc,'k-s','linewidth',1)
plot(SNRout,Pd_H_mc,'k-^','linewidth',1)

h_leg = legend('NMF with correct covariance','NMF with SCM','NMF with NSCM',...
'NMF with KL','NMF with CC','NMF with H');

xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
grid on
str = [str_train,'_',num2str(n),'N','_s',num2str(sigma_t),'.mat'];
save (str); 

