clc
clear
close all
load g_CCIter_KGLRT_1N_s0.5.mat
figure(1);
hold on
plot(SNRout,Pd_SCM_mc,'r','linewidth',2,'markersize',10)
plot(SNRout,Pd_CC_mc,'b','linewidth',2,'markersize',10)
plot(SNRout,Pd_CCIter_mc,'g','linewidth',2,'markersize',10)
plot(SNRout,Pd_ML_mc,'k','linewidth',2,'markersize',10)
load g_CCIter_KGLRT_2N_s0.5.mat
plot(SNRout,Pd_SCM_mc,'r-*','linewidth',2,'markersize',10)
plot(SNRout,Pd_CC_mc,'b-*','linewidth',2,'markersize',10)
plot(SNRout,Pd_CCIter_mc,'g-*','linewidth',2,'markersize',10)
plot(SNRout,Pd_ML_mc,'k-*','linewidth',2,'markersize',10)
h_leg = legend('KGLRT with SCM, K=N',...
'KGLRT with CC, K=N','KGLRT with KA-CE, K=N','KGLRT with ML, K=N',...
'KGLRT with SCM, K=2N',...
'KGLRT with CC, K=2N','KGLRT with KA-CE, K=2N','KGLRT with ML, K=2N');

% h_leg = legend('NMF with correct covariance','NMF with SCM','NMF with NSCM',...
% 'NMF with CCIter','NMF with ML','NMF with CC','NMF with H');

xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
load g_CCIter_AMF_1N_s0.1.mat
figure(2);
hold on
plot(SNRout,Pd_SCM_mc,'r','linewidth',2)
plot(SNRout,Pd_CC_mc,'b','linewidth',2)
plot(SNRout,Pd_CCIter_mc,'g','linewidth',2)
plot(SNRout,Pd_ML_mc,'k','linewidth',2)
load g_CCIter_AMF_2N_s0.1.mat
plot(SNRout,Pd_SCM_mc,'r-*','linewidth',2)
plot(SNRout,Pd_CC_mc,'b-*','linewidth',2)
plot(SNRout,Pd_CCIter_mc,'g-*','linewidth',2)
plot(SNRout,Pd_ML_mc,'k-*','linewidth',2)
h_leg = legend('KGLRT with SCM, K=N',...
'KGLRT with CC, K=N','KGLRT with KA-CE, K=N','KGLRT with ML, K=N',...
'KGLRT with SCM, K=2N',...
'KGLRT with CC, K=2N','KGLRT with KA-CE, K=2N','KGLRT with ML, K=2N');

% h_leg = legend('NMF with correct covariance','NMF with SCM','NMF with NSCM',...
% 'NMF with CCIter','NMF with ML','NMF with CC','NMF with H');

xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
grid on
box on