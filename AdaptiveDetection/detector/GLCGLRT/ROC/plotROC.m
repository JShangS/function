clc
clear 
close all
load Pd_CLGLRT2_ROC22Kmu1lambda3s0.1o1_p.mat
labeltsize=20;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=2;
mkft = 10;

figure
hold on
%%5dB
plot(PFA,Pd_CLGLRT_Mlti_mc(:,1),'k-s','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRTCC_Mlti_mc(:,1),'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRT_Mlti_mc(:,1),'b->','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRTNSCM_Mlti_mc(:,1),'c-*','linewidth',linewide1,'MarkerSize',mkft)
%%10dB
plot(PFA,Pd_CLGLRT_Mlti_mc(:,2),'k-s','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRTCC_Mlti_mc(:,2),'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRT_Mlti_mc(:,2),'b->','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRTNSCM_Mlti_mc(:,2),'c-*','linewidth',linewide1,'MarkerSize',mkft)
%%15dB
plot(PFA,Pd_CLGLRT_Mlti_mc(:,3),'k-s','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRTCC_Mlti_mc(:,3),'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRT_Mlti_mc(:,3),'b->','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRTNSCM_Mlti_mc(:,3),'c-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('GLC-GLRT,K=2N','1S-GLRT with CC,K=2N','1S-GLRT with SCM,K=2N',...
               '1S-GLRT with NSCM,K=2N');
xlabel('P_f_a','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
set(gca,'Xscale','log')
axis([1e-4,1e-1,0,1])
grid on
box on