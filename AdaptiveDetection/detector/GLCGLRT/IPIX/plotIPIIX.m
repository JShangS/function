clc
clear 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Êµ²âÊý¾ÝIPIX%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Pd_CLGLRT_IPIX4_1Ks0.mat
figure
hold on
SNRout = SNRout -15;
plot(SNRout,Pd_CLGLRT_mc,'k-o','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTCC_mc,'g-o','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRT_mc,'b-o','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-o','linewidth',2,'markersize',10);
load Pd_CLGLRT_IPIX4_2Ks0.mat
SNRout = SNRout -15;
plot(SNRout,Pd_CLGLRT_mc,'k-s','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTCC_mc,'g-s','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRT_mc,'b-s','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-s','linewidth',2,'markersize',10);
load Pd_CLGLRT_IPIX4_3Ks0.mat
SNRout = SNRout -15;
plot(SNRout,Pd_CLGLRT_mc,'k-*','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTCC_mc,'g-*','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRT_mc,'b-*','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-*','linewidth',2,'markersize',10);
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(gcf,'Position',[700 0 1200 1000])
grid on
box on
h_leg = legend('GLC-GLRT£¬K=8','1S-GLRT with CC£¬K=8',...
               '1S-GLRT with SCM£¬K=8','1S-GLRT with NSCM£¬K=8',...
                'GLC-GLRT£¬K=16','1S-GLRT with CC£¬K=16',...
                '1S-GLRT with SCM£¬K=16','1S-GLRT with NSCM£¬K=16',...
                'GLC-GLRT£¬K=24','1S-GLRT with CC£¬K=24',...
                '1S-GLRT with SCM£¬K=24','1S-GLRT with NSCM£¬K=24');
set(h_leg,'Location','SouthEast')