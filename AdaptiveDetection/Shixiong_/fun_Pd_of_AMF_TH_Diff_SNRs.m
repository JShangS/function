function fun_Pd_of_AMF_TH_Diff_SNRs(N,L,SNRout,PFA,cos2)
% clear all;  close all; clc
% 刘维建
% 2012.11.21.（2.12.05.16）
% AMF的检测概率和虚警概率计算

%%%%%%%%%%%%%%% 可 变 参 数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=12; % N=4;
L=round(2*N);        % 训练样本数
cos2=0.9;     % 广义夹角余弦
PFA=1E-6;
SNRout=0:1:20;
K=L-N+2;
sin2=1-cos2;  % 广义夹角正弦
SNRnum=10.^(SNRout/10);
Beta=linspace(0,1,200);
step_beta=Beta(2)-Beta(1);
PD_Beta=factorial(L)/factorial(K-1)/factorial(N-2)*Beta.^(K-1).*(1-Beta).^(N-2);
%% 检测门限的计算
eta_amf=linspace(10*eps,24,1e4);
for i=1:length(eta_amf)
    Temp=zeros(size(Beta));
    for m=0:K-2
        Temp=Temp+nchoosek(K-1,m+1)*(Beta*eta_amf(i)).^m;
    end
    PFA_amf_beta=1-(Beta*eta_amf(i))./(1+Beta*eta_amf(i)).^(K-1).*Temp;
    PFA_amf(i)=sum(PFA_amf_beta.*PD_Beta)*step_beta;
end
disp('the threshold is determined')
Diff=abs(PFA_amf-PFA);
[Min, Index]=min(Diff);
figure;plot(eta_amf,log10(PFA_amf),'k.-');% axis([0,Max_TH,0 5*PFA])
hold on; plot([min(eta_amf),max(eta_amf)],[log10(PFA),log10(PFA)],'m--','linewidth',2)
xlabel('门限\eta的区间');ylabel('log10(PFA)');legend('不同门限对应的PFA','设定的PFA')
if Min>PFA*0.1
   error('the computation of the threshold is wrong'); 
else
    TH_AMF_TH=eta_amf(Index); %%%%%%%%%  检测门限  %%%%%%%%%% L=12,Th=1.0397; L=24,Th=1.1421
end
%% PD的计算
for i=1:length(SNRout)
    deta_beta2=Beta*SNRnum(i)*cos2; % delt_Beta2=Beta*SNRnum(el)*cos2;
    c=SNRnum(i)*sin2;
    PDF_Beta=fun_PDF_Noncentral_Beta_FiniteSum(L-N+2,N-1,c,Beta); 
    temp_amf=zeros(size(Beta));
    for m=0:K-2
        temp_amf=temp_amf+nchoosek(K-1,m+1)*(TH_AMF_TH*Beta).^m.*fun_IG(m+1,deta_beta2./(1+Beta*TH_AMF_TH));
    end
    PD_amf_beta=1-(TH_AMF_TH*Beta)./(1+Beta*TH_AMF_TH).^(K-1).*temp_amf;    
    Pd_AMF(i)=sum(PD_amf_beta.*PDF_Beta)*step_beta;
    i
end
figure;plot(SNRout,Pd_AMF,'ro'); hold on;grid on

clear eta_amf Diff PFA_amf temp_amf weight cos2_tmpt UU 
cd  D:\MATLAB\Mat数据\Monograph\Chp03 
tmpt=L/N;
if  mod(tmpt,2) % L不是N的整数倍
    eval(['save Pd_AMF_TH_Diff_SNRs_N' num2str(N) '_L'  num2str(L)   '_SNR'  num2str(min(SNRout)) 't'   num2str(max(SNRout)) '_cos2e'   num2str(cos2) '_PFAm' num2str(-log10(PFA)) '.mat'])
else
   eval(['save Pd_AMF_TH_Diff_SNRs_N' num2str(N) '_L' num2str(tmpt) 'N_SNR'  num2str(min(SNRout)) 't'   num2str(max(SNRout)) '_cos2e'   num2str(cos2) '_PFAm' num2str(-log10(PFA)) '.mat'])
end

if  mod(tmpt,2) % L不是N的整数倍
    eval(['load Pd_8RankOneDetectors_MC_Diff_SNRs_N' num2str(N) '_L'  num2str(L)   '_SNR'  num2str(min(SNRout)) 't'   num2str(max(SNRout)) '_cos2e'   num2str(cos2) '_PFAm' num2str(-log10(PFA)) '.mat'])
else
   eval(['load Pd_8RankOneDetectors_MC_Diff_SNRs_N' num2str(N) '_L' num2str(tmpt) 'N_SNR'  num2str(min(SNRout)) 't'   num2str(max(SNRout)) '_cos2e'   num2str(cos2) '_PFAm' num2str(-log10(PFA)) '.mat'])
end
plot(SNRout,Pd_AMF_mc,'b--')

