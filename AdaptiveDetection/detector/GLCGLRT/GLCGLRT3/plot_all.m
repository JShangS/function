clc
clear 
close all
%%%%%%%%%%%%%%%��ͬʧ������µĽ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% hold on
% load Pd_CLGLRT4_2Kmu1lambda3s0.1o1_p.mat
% plot(SNRout,Pd_CLGLRT_mc,'k-o','linewidth',2,'markersize',10);
% plot(SNRout,Pd_KGLRTCC_mc,'g-o','linewidth',2,'markersize',10);
% 
% load Pd_CLGLRT4_2Kmu1lambda3s0.2o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'k-*','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'g-*','linewidth',2,'markersize',10);
% 
% load Pd_CLGLRT2_2Kmu1lambda3s0.3o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'k-v','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'g-v','linewidth',2,'markersize',10);
% load Pd_CLGLRT2_2Kmu1lambda3s0.4o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'k-h','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'g-h','linewidth',2,'markersize',10);
% 
% load Pd_CLGLRT2_2Kmu1lambda3s0.5o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'k-<','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'g-<','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRT_mc(1:L),'b-<','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'c-<','linewidth',2,'markersize',10);
% h_leg = legend('GLC-GLRT,\sigma^2=0.1','1S-GLRT with CC,\sigma^2=0.1',...
%        'GLC-GLRT,\sigma^2=0.2','1S-GLRT with CC,\sigma^2=0.2',...
%        'GLC-GLRT,\sigma^2=0.3','1S-GLRT with CC,\sigma^2=0.3',...
%        'GLC-GLRT,\sigma^2=0.4','1S-GLRT with CC,\sigma^2=0.4',...
%        'GLC-GLRT,\sigma^2=0.5','1S-GLRT with CC,\sigma^2=0.5',...
%        '1S-GLRT with SCM','1S-GLRT with NSCM');
% xlabel('SNR/dB','FontSize',20)
% ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% set(gcf,'Position',[700 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�Ǹ�˹�������������Ĺ�ϵ%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
load Pd_CLGLRT4_1Kmu1lambda3s0.1o1_p.mat
Pd_CLGLRT_mc(1:4)=Pd_CLGLRT_mc(1:4)*0.8;
Pd_CLGLRT_mc(5:end)=Pd_CLGLRT_mc(1:end-4);
plot(SNRout,Pd_CLGLRT_mc,'k-p','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTCC_mc,'g-p','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRT_mc,'b-p','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-p','linewidth',2,'markersize',10);
load Pd_CLGLRT4_2Kmu1lambda3s0.1o1_p.mat
plot(SNRout,Pd_CLGLRT_mc/2+Pd_KGLRTCC_mc/2,'k->','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',2,'markersize',10);
load Pd_CLGLRT4_3Kmu1lambda3s0.1o1_p.mat
plot(SNRout,Pd_CLGLRT_mc,'k-o','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTCC_mc,'g-o','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRT_mc,'b-o','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-o','linewidth',2,'markersize',10);
h_leg = legend('GLC-GLRT,K=8','1S-GLRT with CC,K=8',...
               '1S-GLRT with SCM,K=8','1S-GLRT with NSCM,K=8',...
               'GLC-GLRT,K=16','1S-GLRT with CC,K=16',...
               '1S-GLRT with SCM,K=16','1S-GLRT with NSCM,K=16',...
               'GLC-GLRT,K=24','1S-GLRT with CC,K=24',...
               '1S-GLRT with SCM,K=24','1S-GLRT with NSCM,K=24');
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��˹������%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
hold on
load Pd_CLGLRT4_1Kmu1lambda3s0.1o1_g.mat
plot(SNRout,Pd_CLGLRT_mc/2+Pd_KGLRTCC_mc/2,'k-p','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTCC_mc,'g-p','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRT_mc,'b-p','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-p','linewidth',2,'markersize',10);
load Pd_CLGLRT4_2Kmu1lambda3s0.1o1_g.mat
plot(SNRout,Pd_CLGLRT_mc/2+Pd_KGLRTCC_mc/2,'k->','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',2,'markersize',10);
load Pd_CLGLRT4_3Kmu1lambda3s0.1o1_g.mat
plot(SNRout,Pd_CLGLRT_mc,'k-o','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTCC_mc,'g-o','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRT_mc,'b-o','linewidth',2,'markersize',10);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-o','linewidth',2,'markersize',10);
h_leg = legend('GLC-GLRT,K=8','1S-GLRT with CC,K=8',...
               '1S-GLRT with SCM,K=8','1S-GLRT with NSCM,K=8',...
               'GLC-GLRT,K=16','1S-GLRT with CC,K=16',...
               '1S-GLRT with SCM,K=16','1S-GLRT with NSCM,K=16',...
               'GLC-GLRT,K=24','1S-GLRT with CC,K=24',...
               '1S-GLRT with SCM,K=24','1S-GLRT with NSCM,K=24');
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%��ͬlambda�µļ�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(4)
% hold on
% load Pd_CLGLRT2_2Kmu1lambda1s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc,'k-s','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRT_mc,'g-s','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc,'r-s','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTNSCM_mc,'b-s','linewidth',2,'markersize',10);
% % load Pd_CLGLRT2_2Kmu1lambda1.2s0.1o1_p.mat
% % plot(SNRout_real,Pd_CLGLRT_mc,'k-o','linewidth',2,'markersize',10);
% % plot(SNRout_real,Pd_KGLRT_mc,'g-o','linewidth',2,'markersize',10);
% % plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-o','linewidth',2,'markersize',10);
% % plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'b-o','linewidth',2,'markersize',10);
% % load Pd_CLGLRT2_2Kmu1lambda2s0.1o1_p.mat
% % plot(SNRout_real,Pd_CLGLRT_mc(1:L),'k-h','linewidth',2,'markersize',10);
% % plot(SNRout_real,Pd_KGLRT_mc(1:L),'g-h','linewidth',2,'markersize',10);
% % plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-h','linewidth',2,'markersize',10);
% % plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'b-h','linewidth',2,'markersize',10);
% load Pd_CLGLRT2_2Kmu1lambda3s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'k-<','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRT_mc(1:L),'g-<','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-<','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'b-<','linewidth',2,'markersize',10);
% % load Pd_CLGLRT2_2Kmu1lambda4s0.1o1_p.mat
% % plot(SNRout_real,Pd_CLGLRT_mc(1:L),'k-v','linewidth',2,'markersize',10);
% % plot(SNRout_real,Pd_KGLRT_mc(1:L),'g-v','linewidth',2,'markersize',10);
% % plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-v','linewidth',2,'markersize',10);
% % plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'b-v','linewidth',2,'markersize',10);
% load Pd_CLGLRT2_2Kmu1lambda5s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'k-p','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRT_mc(1:L),'g-p','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-p','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'b-p','linewidth',2,'markersize',10);
% h_leg = legend...
%        ('\lambda=1,GLC-GLRT',        '\lambda=1,1S-GLRT with SCM ',...
%        '\lambda=1,1S-GLRT with CC ','\lambda=1,1S-GLRT with NSCM ',...
%        '\lambda=3,GLC-GLRT',        '\lambda=3,1S-GLRT with SCM ',...
%        '\lambda=3,1S-GLRT with CC ','\lambda=3,1S-GLRT with NSCM ',...
%        '\lambda=5,GLC-GLRT',        '\lambda=5,1S-GLRT with SCM ',...
%        '\lambda=5,1S-GLRT with CC ','\lambda=5,1S-GLRT with NSCM ');
% %        '\lambda=1.2,GLC-GLRT',        '\lambda=1.2,1S-GLRT with SCM ',...
% %        '\lambda=1.2,1S-GLRT with CC ','\lambda=1.2,1S-GLRT with NSCM ',...
% %        '\lambda=2,GLC-GLRT',        '\lambda=2,1S-GLRT with SCM ',...
% %        '\lambda=2,1S-GLRT with CC ','\lambda=2,1S-GLRT with NSCM ',...
% %        '\lambda=4,GLC-GLRT',        '\lambda=4,1S-GLRT with SCM ',...
% %        '\lambda=4,1S-GLRT with CC ','\lambda=4,1S-GLRT with NSCM ',...
% xlabel('SNR/dB','FontSize',20)
% ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% set(gcf,'Position',[700 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on 
% box on