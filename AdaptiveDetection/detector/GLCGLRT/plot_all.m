clc
clear 
close all
%%%%%%%%%%%%%%%不同失配情况下的结果%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% hold on
% load Pd_CLGLRT2_2Kmu1lambda3s0.1o1_p.mat
% SNRout_real=0:1:20; % 输出SNR
% L = length(SNRout_real);
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'k-o','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'g-o','linewidth',2,'markersize',10);
% 
% load Pd_CLGLRT2_2Kmu1lambda3s0.2o1_p.mat
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%非高斯环境下数据量的关系%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
load Pd_CLGLRT3_1Kmu1lambda3s0.1o1_p.mat
SNRout_real=0:1:20; % 输出SNR
L = length(SNRout_real);
plot(SNRout_real,Pd_CLGLRT_mc(1:L),'k-p','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'g-p','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRT_mc(1:L),'b-p','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'c-p','linewidth',2,'markersize',10);
% load Pd_CLGLRT2_1.5Kmu1lambda3s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-s','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-s','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRT_mc(1:L),'g-s','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'c-s','linewidth',2,'markersize',10);
load Pd_CLGLRT3_2Kmu1lambda3s0.1o1_p.mat
plot(SNRout_real,Pd_CLGLRT_mc(1:L),'k->','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'g->','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRT_mc(1:L),'b->','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'c->','linewidth',2,'markersize',10);
load Pd_CLGLRT3_3Kmu1lambda3s0.1o1_p.mat
plot(SNRout_real,Pd_CLGLRT_mc(1:L),'k-o','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'g-o','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRT_mc(1:L),'b-o','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'c-o','linewidth',2,'markersize',10);
h_leg = legend('GLC-GLRT，K=16','1S-GLRT with CC，K=16',...
       '1S-GLRT with SCM，K=16','1S-GLRT with NSCM，K=16',...
       'GLC-GLRT，K=32','1S-GLRT with CC，K=32',...
       '1S-GLRT with SCM，K=32','1S-GLRT with NSCM，K=32',...
       'GLC-GLRT，K=48','1S-GLRT with CC，K=48',...
       '1S-GLRT with SCM，K=48','1S-GLRT with NSCM，K=48');
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%高斯环境下%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
hold on
load Pd_CLGLRT2_1Kmu1lambda3s0.1o1_g.mat
SNRout_real=0:1:20; % 输出SNR
L = length(SNRout_real);
plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-s','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-s','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRT_mc(1:L),'g-s','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'k-s','linewidth',2,'markersize',10);
load Pd_CLGLRT2_2Kmu1lambda3s0.1o1_g.mat
plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b->','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r->','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRT_mc(1:L),'g->','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'k->','linewidth',2,'markersize',10);
load Pd_CLGLRT2_4Kmu1lambda3s0.1o1_g.mat
plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-o','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-o','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRT_mc(1:L),'g-o','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'k-o','linewidth',2,'markersize',10);
h_leg = legend('GLC-GLRT,K=16','1S-GLRT with CC,K=16',...
               '1S-GLRT with SCM,K=16','1S-GLRT with NSCM,K=16',...
               'GLC-GLRT,K=32','1S-GLRT with CC,K=32',...
               '1S-GLRT with SCM,K=32','1S-GLRT with NSCM,K=32',...
               'GLC-GLRT,K=64','1S-GLRT with CC,K=64',...
               '1S-GLRT with SCM,K=64','1S-GLRT with NSCM,K=64');
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%不同lambda下的检测概率%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(4)
% hold on
% load Pd_CLGLRT2_2Kmu1lambda1s0.1o1_p.mat
% SNRout_real=0:1:20; % 输出SNR
% L = length(SNRout_real);
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'k-s','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRT_mc(1:L),'g-s','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-s','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'b-s','linewidth',2,'markersize',10);
% % load Pd_CLGLRT2_2Kmu1lambda1.2s0.1o1_p.mat
% % plot(SNRout_real,Pd_CLGLRT_mc(1:L),'k-o','linewidth',2,'markersize',10);
% % plot(SNRout_real,Pd_KGLRT_mc(1:L),'g-o','linewidth',2,'markersize',10);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%带有AML比对的失配效果%%%%%%%%%%%%%%%%%%%%
% figure(5)
% hold on
% SNRout_real=0:1:20; % 输出SNR
% L = length(SNRout_real);
% load Pd_CLGLRT2_2Kmu1lambda3s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-o','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-o','linewidth',2,'markersize',10);
% % 
% load Pd_CLGLRT_2Kmu1lambda4s0.2o1_pML.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-*','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-*','linewidth',2,'markersize',10);
% 
% load Pd_CLGLRT_2Kmu1lambda4s0.3o1_pML.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-v','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-v','linewidth',2,'markersize',10);
% 
% load Pd_CLGLRT_2Kmu1lambda4s0.4o1_pML.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-h','linewidth',2,'markersize',10);
% % plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-h','linewidth',2,'markersize',10);
% % 
% load Pd_CLGLRT2_2Kmu1lambda3s0.5o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-<','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-<','linewidth',2,'markersize',10);
% % 
% plot(SNRout_real,Pd_KGLRT_mc(1:L),'g-<','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_AMF_mc(1:L),'g-o','linewidth',2,'markersize',10);
% 
% legend('CLGLRT,\sigma^2=0.1','GLRTCC,\sigma^2=0.1',...
%        'CLGLRT,\sigma^2=0.5','GLRTCC,\sigma^2=0.5')
% %    'CLGLRT,\sigma^2=0.2','GLRTCC,\sigma^2=0.2','GLRTML,\sigma^2=0.2','AMFCC,\sigma^2=0.2','AMFML,\sigma^2=0.2',...
% %        'CLGLRT,\sigma^2=0.3','GLRTCC,\sigma^2=0.3','GLRTML,\sigma^2=0.3','AMFCC,\sigma^2=0.3','AMFML,\sigma^2=0.3',...
% %        'CLGLRT,\sigma^2=0.4','GLRTCC,\sigma^2=0.4','GLRTML,\sigma^2=0.4','AMFCC,\sigma^2=0.4','AMFML,\sigma^2=0.4',...
% xlabel('SNR/dB','FontSize',20)
% ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% set(gcf,'Position',[600 0 900 800])
% grid on
% box on
% figure(1)
% hold on
% load Pd_CLGLRT2_1Kmu1lambda2s0.1o1_p.mat
% SNRout_real=0:1:20; % 输出SNR
% L = length(SNRout_real);
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-<','linewidth',2,'markersize',10);
% % load Pd_CLGLRT2_1.5Kmu1lambda4s0.1o1_p.mat
% % plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-s','linewidth',2,'markersize',10);
% load Pd_CLGLRT2_2Kmu1lambda2s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-o','linewidth',2,'markersize',10);
% load Pd_CLGLRT2_4Kmu1lambda2s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-p','linewidth',2,'markersize',10);
% % legend('1','1.5','2','4')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%实测数据PhaseOne%%%%%%%%%%%%%%%%%%%%%%
% load Pd_CLGLRT2_real1Ks0.1.mat
% figure
% hold on
% plot(SNRout-15,Pd_CLGLRT_mc,'k-o','linewidth',2,'markersize',10);
% plot(SNRout-15,Pd_KGLRTCC_mc,'g-o','linewidth',2,'markersize',10);
% plot(SNRout-15,Pd_KGLRT_mc,'b-o','linewidth',2,'markersize',10);
% plot(SNRout-15,Pd_KGLRTNSCM_mc,'c-o','linewidth',2,'markersize',10);
% load Pd_CLGLRT2_real2Ks0.1.mat
% plot(SNRout-10,Pd_CLGLRT_mc,'k-s','linewidth',2,'markersize',10);
% plot(SNRout-10,Pd_KGLRTCC_mc,'g-s','linewidth',2,'markersize',10);
% plot(SNRout-10,Pd_KGLRT_mc,'b-s','linewidth',2,'markersize',10);
% plot(SNRout-10,Pd_KGLRTNSCM_mc,'c-s','linewidth',2,'markersize',10);
% load Pd_CLGLRT2_real3Ks0.1.mat
% plot(SNRout-10,Pd_CLGLRT_mc,'k-*','linewidth',2,'markersize',10);
% plot(SNRout-10,Pd_KGLRTCC_mc,'g-*','linewidth',2,'markersize',10);
% plot(SNRout-10,Pd_KGLRT_mc,'b-*','linewidth',2,'markersize',10);
% plot(SNRout-10,Pd_KGLRTNSCM_mc,'c-*','linewidth',2,'markersize',10);
% xlabel('SNR/dB','FontSize',20)
% ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% set(gcf,'Position',[700 0 1200 1000])
% grid on
% box on
% h_leg = legend('GLC-GLRT，K=16','1S-GLRT with CC，K=16',...
%                '1S-GLRT with SCM，K=16','1S-GLRT with NSCM，K=16',...
%                 'GLC-GLRT，K=32','1S-GLRT with CC，K=32',...
%                 '1S-GLRT with SCM，K=32','1S-GLRT with NSCM，K=32',...
%                 'GLC-GLRT，K=48','1S-GLRT with CC，K=48',...
%                 '1S-GLRT with SCM，K=48','1S-GLRT with NSCM，K=48')
% set(h_leg,'Location','SouthEast')