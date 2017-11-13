%%%ʵ��һ������ɫ���ص�GLRT�����
clc
clear 
close all
Na = 4;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
SNRout=0:1:20; % ���SNR
cos2=0.9;
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
rou = 0.95;  %%Э����������ɵĳ�������
rouR = zeros(N,N);  %%��ʵ���Ӳ�Э����
L=round(1*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% ϵͳ����ʸ��
for i=1:N
    for j=1:N
        rouR(i,j)=rou^abs(i-j);%*exp(1j*2*pi*abs(i-j)*theta_sig)
    end
end
irouR=inv(rouR);
rouR_abs=abs(rouR);
t = normrnd(1,0.1,N,1);%%0~0.5%%ʧ������
R_KA = rouR.*(t*t');
% R_KA_inv = inv(R_KA);
rouR_half=rouR^0.5;


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
% figure; plot(weight,cos2_tmpt);


% Train = fun_TrainData(N,L,rouR);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
% x0 = fun_TrainData(N,1,rouR); % �����źŽ������Ӳ�������
% x0=s_real+x0;%+pp;    %%%%%%%  ��Ҫ  %%%%%%%%%%%%%



lamda = 1;
mu = 1;
%%%���޼���
h = waitbar(0,'Please wait...');
for i = 1:MonteCarloPfa
    waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
    Train = fun_TrainData(N,L,rouR);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    x0 = fun_TrainData(N,1,rouR); % �����źŽ������Ӳ�������
    %%%%Э�������
    R_SCM = abs(fun_SCM(Train));
    iR_SCM=inv(R_SCM);
    R_ICL1 = (fun_ICL(s,x0,inv(R_KA),inv(R_SCM),lamda,mu,1));
    R_ICL0 = (fun_ICL(s,x0,inv(R_KA),inv(R_SCM),lamda,mu,0));
    iR_ICL1 = inv(R_ICL1);
    iR_ICL0 = inv(R_ICL0);
    %%%�����
    Tamf(i) = abs(s'*iR_SCM*x0)^2/abs(s'*iR_SCM*s);     %%%%%% AMF����wald
    tmp=abs(x0'*iR_SCM*x0);
    Tglrt(i) = Tamf(i)/(1+tmp);                   %%%%%% KGLRT
    a = (s'*iR_ICL1*x0)/(s'*iR_ICL1*s);
    tmp1 = ((x0 - a*s)'*iR_ICL1*(x0 - a*s)+1/mu)^(-lamda-N);
    tmp2 = (x0'*iR_ICL1*x0+1/mu)^(-lamda-N);
    Tclglrt(i) =  abs(tmp1/tmp2);%%%%%% ɫ���ص�GLRT
end
close(h)
TGLRT=sort(Tglrt,'descend');
TCLGLRT=sort(Tclglrt,'descend');
Th_GLRT=(TGLRT(floor(MonteCarloPfa*PFA-1))+TGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_CLGLRT=(TCLGLRT(floor(MonteCarloPfa*PFA-1))+TCLGLRT(floor(MonteCarloPfa*PFA)))/2;

counter_glrt=0;
counter_clglrt=0;
Pd_GLRT_mc = zeros(1,length(SNRout));
Pd_CLGLRT_mc = zeros(1,length(SNRout));
alpha=sqrt(SNRnum/abs(s_real'*irouR*s_real)); % ����SNR=|alpha|^2*s'*R^(-1)*s���|alpha|
h = waitbar(0,'Please wait...');
for m=1:length(SNRout)
    
    for i=1:MonteCarloPd 
        waitbar(((m-1)*MonteCarloPd+i)/length(SNRout)/MonteCarloPd,h,sprintf([num2str(((m-1)*MonteCarloPd+i)/length(SNRout)/MonteCarloPd*100),'%%']));
        Train = fun_TrainData(N,L,rouR);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
        x0 = fun_TrainData(N,1,rouR); % �����źŽ������Ӳ�������
        %%%%Э�������
        R_SCM = abs(fun_SCM(Train));
        iR_SCM=inv(R_SCM);
        x0=alpha(m)*s_real+x0;%+pp;    %%%%%%%  ��Ҫ  %%%%%%%%%%%%%
        R_ICL1 = (fun_ICL(s,x0,inv(R_KA),inv(R_SCM),lamda,mu,1));
        R_ICL0 = (fun_ICL(s,x0,inv(R_KA),inv(R_SCM),lamda,mu,0));
        iR_ICL1 = inv(R_ICL1);
        iR_ICL0 = inv(R_ICL0);
        %%%�����
        Tamf = abs(s'*iR_SCM*x0)^2/abs(s'*iR_SCM*s);     %%%%%% AMF����wald
        tmp=abs(x0'*iR_SCM*x0);
        Tglrt = Tamf/(1+tmp);                   %%%%%% KGLRT
%         a = (s'*iR_ICL1*x0)/(s'*iR_ICL1*s);
%         tmp1 = ((x0 - a*s)'*iR_ICL1*(x0 - a*s)+1/mu)^(-lamda-N);
%         tmp2 = (x0'*iR_ICL1*x0+1/mu)^(-lamda-N);
%         Tclglrt =  abs(tmp1/tmp2);%%%%%% ɫ���ص�GLRT
        if Tglrt>Th_GLRT;         counter_glrt=counter_glrt+1;          end            
%         if Tclglrt>Th_CLGLRT;      counter_clglrt=counter_clglrt+1;        end   
    end
    Pd_GLRT_mc(m)=counter_glrt/MonteCarloPd;          counter_glrt=0;
%     Pd_CLGLRT_mc(m)=counter_clglrt/MonteCarloPd;        counter_clglrt=0;
end
close(h)
figure(2);
hold on
plot(SNRout,Pd_GLRT_mc,'b-+','linewidth',2)
plot(SNRout,Pd_CLGLRT_mc,'g.-')
