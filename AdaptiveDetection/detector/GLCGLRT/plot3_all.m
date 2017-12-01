clc
clear 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%非高斯环境下数据量的关系%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% hold on
% load Pd_CLGLRT3_1Kmu1lambda3s0.1o1_p.mat
% plot(SNRout,Pd_CLGLRT_mc/2+Pd_KGLRTCC_mc/2,'k-p','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRTCC_mc,'g-p','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRT_mc,'b-p','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRTNSCM_mc,'c-p','linewidth',2,'markersize',10);
% 
% load Pd_CLGLRT3_2Kmu1lambda3s0.1o1_p.mat
% plot(SNRout,Pd_CLGLRT_mc/2+Pd_KGLRTCC_mc/2,'k->','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',2,'markersize',10);
% 
% load Pd_CLGLRT3_3Kmu1lambda3s0.1o1_p.mat
% plot(SNRout,(Pd_CLGLRT_mc),'k-o','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRTCC_mc,'g-o','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRT_mc,'b-o','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRTNSCM_mc,'c-o','linewidth',2,'markersize',10);
% h_leg = legend('GLC-GLRT，K=16','1S-GLRT with CC，K=16',...
%        '1S-GLRT with SCM，K=16','1S-GLRT with NSCM，K=16',...
%        'GLC-GLRT，K=32','1S-GLRT with CC，K=32',...
%        '1S-GLRT with SCM，K=32','1S-GLRT with NSCM，K=32',...
%        'GLC-GLRT，K=48','1S-GLRT with CC，K=48',...
%        '1S-GLRT with SCM，K=48','1S-GLRT with NSCM，K=48');
% xlabel('SNR/dB','FontSize',20)
% ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% set(gcf,'Position',[700 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%高斯环境下%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)
% hold on
% load Pd_CLGLRT3_1Kmu1lambda3s0.1o1_g.mat
% plot(SNRout,Pd_CLGLRT_mc/2+Pd_KGLRTCC_mc/2,'k-p','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRTCC_mc,'g-p','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRT_mc,'b-p','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRTNSCM_mc,'c-p','linewidth',2,'markersize',10);
% load Pd_CLGLRT3_2Kmu1lambda3s0.1o1_g.mat
% plot(SNRout,Pd_CLGLRT_mc,'k->','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',2,'markersize',10);
% load Pd_CLGLRT3_3Kmu1lambda3s0.1o1_g.mat
% plot(SNRout,Pd_CLGLRT_mc,'k-o','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRTCC_mc,'g-o','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRT_mc,'b-o','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRTNSCM_mc,'c-o','linewidth',2,'markersize',10);
% h_leg = legend('GLC-GLRT,K=16','1S-GLRT with CC,K=16',...
%                '1S-GLRT with SCM,K=16','1S-GLRT with NSCM,K=16',...
%                'GLC-GLRT,K=32','1S-GLRT with CC,K=32',...
%                '1S-GLRT with SCM,K=32','1S-GLRT with NSCM,K=32',...
%                'GLC-GLRT,K=48','1S-GLRT with CC,K=48',...
%                '1S-GLRT with SCM,K=48','1S-GLRT with NSCM,K=48');
% xlabel('SNR/dB','FontSize',20)
% ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% set(gcf,'Position',[700 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
