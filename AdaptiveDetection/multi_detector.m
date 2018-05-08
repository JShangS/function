clc
clear 
close all
% ��ά����2017.09.15
% ���ּ�����ļ�����ܷ���
%%
Na=2;     % ��Ԫ��
Np=4;     % ������
N=Na*Np;
L=round(2*N); 
SNRout=0:1:25; % ���SNR
cos2=1;%%%ʧ�����
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=round(1/PFA*100);
MonteCarloPd=1e4;
R=zeros(N);
% vt=randn(N,1)+1i*randn(N,1);
theta_sig=0.1;
nn=0:N-1;
vt=exp(-1i*2*pi*nn*theta_sig)'; %%%%%% ϵͳ����ʸ��
rhoR=0.95;
for i=1:N
    for j=1:N
        R(i,j)=rhoR^abs(i-j)*exp(1j*2*pi*(i-j)*0.1);
    end
end
iR=inv(R);
R_half=R^0.5;
% [Ddescend,Index]=sort(diagD,'descend'); figure; plot(10*log10(abs(Ddescend)),'b-x')
[UU,SS,VV]=svd(iR*vt);
vt_v=UU(:,2); %%%%%% ��vt�ڰ׻��ռ�����������vt^H*iR*vt_v==0
weight=linspace(0,1,300);
for i=1:length(weight)
    vt_tmpt=weight(i)*vt+(1-weight(i))*vt_v;
    cos2_tmpt(i)=abs(vt_tmpt'*iR*vt).^2/abs(vt_tmpt'*iR*vt_tmpt*vt'*iR*vt);
end
[Min, Index]=min(abs(cos2-cos2_tmpt));
Weight=weight(Index)
vt_real=Weight*vt+(1-Weight)*vt_v;
figure; plot(weight,cos2_tmpt);
%% ����������
Tamf = zeros(1,MonteCarloPfa);
Tglrt = zeros(1,MonteCarloPfa);
Tace = zeros(1,MonteCarloPfa);
Tabort = zeros(1,MonteCarloPfa);
Twabort = zeros(1,MonteCarloPfa);
Tprao = zeros(1,MonteCarloPfa);
Tdnamf = zeros(1,MonteCarloPfa);
Taed = zeros(1,MonteCarloPfa);
tic    
% h = waitbar(0,'Please wait...');
parfor i=1:MonteCarloPfa
%     waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
%     X=(randn(N,L)+1i*randn(N,L))/sqrt(2); % ��������Ϊ1�ĸ���˹������ % Rwhite1=1/snapshot1*X1*X1'; eig(Rwhite1); % round(mean(abs(eig(Rwhite1)))) == 1
%     S=(R_half*X)*(R_half*X)'; % ��L��ѵ���������Ƶ��Ӳ���������Э�������(Rhalf*X��ʾ���յ�L��ѵ������)
    Train = fun_TrainData('g',N,L,R,3,1,1);
    S = fun_SCMN(Train);
    iS=inv(S);
%     W=(randn(N,1)+1i*randn(N,1))/sqrt(2); % 1i == -i
%     x=R_half*W;%+pp; % noise=(randn(N,1)+j*randn(N,1))/sqrt(2);  % �����źŽ������Ӳ�������
    x = fun_TrainData('g',N,1,R,3,1,1);
    Tamf(i)=abs(vt'*iS*x)^2/abs(vt'*iS*vt);     %%%%%% AMF����wald
    tmp=abs(x'*iS*x);
    Tglrt(i)=Tamf(i)/(1+tmp);                   %%%%%% KGLRT
    Tace(i)=Tamf(i)/tmp;                        %%%%%% ACE
    Tabort(i)=(1+Tamf(i))/(2+tmp);              %%%%%% ABORT  % eq.(16) ���ͳ����
    Twabort(i)=1/(1+tmp)/(1-Tglrt(i))^2;        %%%%%% ABORT  % �����������е�eq.(18)
    Tace_bar=Tace(i)/(1-Tace(i));
    Tprao(i)=Tglrt(i)^2/(Tamf(i)*(1-Tglrt(i))); %%%%%% DMRao
    Tdnamf(i)=Tace_bar/tmp;                     %%%%%% DNAMF  % eq.(24) ���ͳ����
    Taed(i)=tmp;                                %%%%%% ��������� 
    Tmtd_t = fft(x);
    Tmtd(i) = abs(Tmtd_t(2));
end
TAMF=sort(Tamf,'descend');
TACE=sort(Tace,'descend');
TKGLRT=sort(Tglrt,'descend');
TABORT=sort(Tabort,'descend');
TWABORT=sort(Twabort,'descend');
TDMRao=sort(Tprao,'descend');
TDNAMF=sort(Tdnamf,'descend');
TAED=sort(Taed,'descend');
TMTD = sort(Tmtd,'descend');

Th_AMF=(TAMF(floor(MonteCarloPfa*PFA-1))+TAMF(floor(MonteCarloPfa*PFA)))/2;
Th_ACE=(TACE(floor(MonteCarloPfa*PFA-1))+TACE(floor(MonteCarloPfa*PFA)))/2;
Th_KGLRT=(TKGLRT(floor(MonteCarloPfa*PFA-1))+TKGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_ABORT=(TABORT(floor(MonteCarloPfa*PFA-1))+TABORT(floor(MonteCarloPfa*PFA)))/2;
Th_WABORT=(TWABORT(floor(MonteCarloPfa*PFA-1))+TWABORT(floor(MonteCarloPfa*PFA)))/2;
Th_DMRao=(TDMRao(floor(MonteCarloPfa*PFA-1))+TDMRao(floor(MonteCarloPfa*PFA)))/2;
Th_DNAMF=(TDNAMF(floor(MonteCarloPfa*PFA-1))+TDNAMF(floor(MonteCarloPfa*PFA)))/2;
Th_AED=(TAED(floor(MonteCarloPfa*PFA-1))+TAED(floor(MonteCarloPfa*PFA)))/2;
Th_MTD=(TMTD(floor(MonteCarloPfa*PFA-1))+TMTD(floor(MonteCarloPfa*PFA)))/2;
toc
% alpha=sqrt(SNRnum/abs(vt'*invR*vt)); % ����SNR=|alpha|^2*s'*R^(-1)*s���|alpha|
a=0;b=2*pi;
%% ���������
tic
counter_amf=0;
counter_ace=0;
counter_glrt=0;
counter_abort=0;
counter_wabort=0;
counter_prao=0;
counter_dnamf=0;
counter_aed=0;
counter_mtd = 0;
alpha=sqrt(SNRnum/abs(vt'*iR*vt)); % ����SNR=|alpha|^2*s'*R^(-1)*s���|alpha|
Pd_AMF_mc = zeros(1,length(SNRout));
Pd_KGLRT_mc = zeros(1,length(SNRout));
Pd_ACE_mc = zeros(1,length(SNRout));
Pd_ABORT_mc = zeros(1,length(SNRout));
Pd_WABORT_mc = zeros(1,length(SNRout));
Pd_DMRao_mc = zeros(1,length(SNRout));
Pd_DNAMF_mc = zeros(1,length(SNRout));
Pd_AED_mc = zeros(1,length(SNRout));
Pd_MTD_mc = zeros(1,length(SNRout));
h = waitbar(0,'Please wait...');
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
    parfor i=1:MonteCarloPd
%         X=(randn(N,L)+1i*randn(N,L))/sqrt(2); % ��������Ϊ1�ĸ���˹������ % Rwhite1=1/snapshot1*X1*X1'; eig(Rwhite1); % round(mean(abs(eig(Rwhite1)))) == 1
%         S=(R_half*X)*(R_half*X)'; % ��L��ѵ���������Ƶ��Ӳ���������Э�������(Rhalf*X��ʾ���յ�L��ѵ������)
        Train = fun_TrainData('g',N,L,R,3,1,1);
        S = fun_SCMN(Train);
        iS=inv(S);
        W=(randn(N,1)+1i*randn(N,1))/sqrt(2); % 1i == -i
    %     THETA=a+(b-a)*rand; % ����0--2*pi�ľ��ȷֲ������λ
    %         x=alpha(m)*exp(1i*THETA)*vt+Clutter;%+pp;
        x = fun_TrainData('g',N,1,R,3,1,1);
        x=alpha(m)*vt+x;%+pp;    %%%%%%%  ��Ҫ  %%%%%%%%%%%%%
        Tamf=abs(vt'*iS*x)^2/abs(vt'*iS*vt);  %%%%%% AMF
        tmp=abs(x'*iS*x);
        Tglrt=Tamf/(1+tmp);                   %%%%%% KGLRT
        Tace=Tamf/tmp;                        %%%%%% ACE
        Tabort=(1+Tamf)/(2+tmp);              %%%%%% ABORT  % eq.(16) ���ͳ����
        Twabort=1/(1+tmp)/(1-Tglrt)^2;        %%%%%% ABORT  % �����������е�eq.(18)
        Tace_bar=Tace/(1-Tace);
        Tprao=Tglrt^2/(Tamf*(1-Tglrt));       %%%%%% DMRao
        Tdnamf=Tace_bar/tmp;                  %%%%%% DNAMF  % eq.(24) ���ͳ����
        Taed=tmp;                             %%%%%% ���������  
        Tmtd_t = fft(x);
        Tmtd = abs(Tmtd_t(2));
        if Tamf>Th_AMF;         counter_amf=counter_amf+1;          end            
        if Tglrt>Th_KGLRT;      counter_glrt=counter_glrt+1;        end                
        if Tace>Th_ACE;         counter_ace=counter_ace+1;          end          
        if Tabort>Th_ABORT;     counter_abort=counter_abort+1;      end            
        if Twabort>Th_WABORT;   counter_wabort=counter_wabort+1;    end        
        if Tprao>Th_DMRao;      counter_prao=counter_prao+1;        end          
        if Tdnamf>Th_DNAMF;     counter_dnamf=counter_dnamf+1;      end          
        if Taed>Th_AED;         counter_aed=counter_aed+1;          end
        if Tmtd>Th_MTD;         counter_mtd=counter_mtd+1;          end
    end
    Pd_AMF_mc(m)=counter_amf/MonteCarloPd;          counter_amf=0;
    Pd_KGLRT_mc(m)=counter_glrt/MonteCarloPd;        counter_glrt=0;
    Pd_ACE_mc(m)=counter_ace/MonteCarloPd;          counter_ace=0;
    Pd_ABORT_mc(m)=counter_abort/MonteCarloPd;      counter_abort=0;
    Pd_WABORT_mc(m)=counter_wabort/MonteCarloPd;    counter_wabort=0;
    Pd_DMRao_mc(m)=counter_prao/MonteCarloPd;        counter_prao=0;
    Pd_DNAMF_mc(m)=counter_dnamf/MonteCarloPd;      counter_dnamf=0;
    Pd_AED_mc(m)=counter_aed/MonteCarloPd;          counter_aed=0;  
    Pd_MTD_mc(m)=counter_mtd/MonteCarloPd;          counter_mtd=0; 
    m
end
close(h)
toc

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
plot(SNRout,Pd_MTD_mc,'b-s','linewidth',2); 
legend('KGLRT','AMF/DMwald','ACE','ABORT','WABORT','DMRao','AED','DNAMF','MTD')
% legend({'KGLRT','AMF/DMwald','DMRao'},'FontSize',20)
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
grid on
% axis([5 25 0 1])
% clear TAMF TACE TKGLRT TABORT TWABORT TDMRao TDNAMF   Tamf Tace Tglrt Tabort Twabort Tprao Tdnamf Taed TAED X
% cd  D:\MATLAB\Mat����\Monograph\Chp03 
% tmpt=L/N;
% if  mod(tmpt,2) % L����N��������
%     eval(['save Pd_8RankOneDetectors_MC_Diff_SNRs_N' num2str(N) '_L'  num2str(L)   '_SNR'  num2str(min(SNRout)) 't'   num2str(max(SNRout)) '_cos2e'   num2str(cos2)  '_PFAm' num2str(-log10(PFA)) '.mat'])
% else
%    eval(['save Pd_8RankOneDetectors_MC_Diff_SNRs_N' num2str(N) '_L' num2str(tmpt) 'N_SNR'  num2str(min(SNRout)) 't'   num2str(max(SNRout)) '_cos2e'   num2str(cos2) '_PFAm' num2str(-log10(PFA)) '.mat'])
% end