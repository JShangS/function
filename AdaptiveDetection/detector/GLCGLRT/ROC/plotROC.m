clc
clear 
close all
load Pd_CLGLRT2_ROC22Kmu1lambda3s0.1o1_p.mat
figure
hold on
%%5dB
plot(PFA,Pd_CLGLRT_Mlti_mc(:,1),'k-s','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRTCC_Mlti_mc(:,1),'g-o','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRT_Mlti_mc(:,1),'b->','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRTNSCM_Mlti_mc(:,1),'c-*','linewidth',2,'MarkerSize',15)
%%10dB
plot(PFA,Pd_CLGLRT_Mlti_mc(:,2),'k-s','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRTCC_Mlti_mc(:,2),'g-o','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRT_Mlti_mc(:,2),'b->','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRTNSCM_Mlti_mc(:,2),'c-*','linewidth',2,'MarkerSize',15)
%%15dB
plot(PFA,Pd_CLGLRT_Mlti_mc(:,3),'k-s','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRTCC_Mlti_mc(:,3),'g-o','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRT_Mlti_mc(:,3),'b->','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_KGLRTNSCM_Mlti_mc(:,3),'c-*','linewidth',2,'MarkerSize',15)
h_leg = legend('GLC-GLRT,K=16','1S-GLRT with CC,K=16','1S-GLRT with SCM,K=16',...
               '1S-GLRT with NSCM,K=16');
xlabel('Pfa','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
set(gcf,'Position',[700 0 1200 1000])
grid on
box on