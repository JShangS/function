clc
clear 
close all
labeltsize=20;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=2;
mkft = 10;
%%%%%%%%%%%%%%%不同失配情况下的结果%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
load Pd_CLGLRT4_2Kmu1lambda3s0.1o1_p.mat
plot(SNRout,Pd_CLGLRT_mc,'k-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-o','linewidth',linewide1,'markersize',mkft);

% load Pd_CLGLRT4_2Kmu1lambda3s0.2o1_p.mat
% plot(SNRout,Pd_CLGLRT_mc,'k-*','linewidth',linewide1,'markersize',mkft);
% plot(SNRout,Pd_KGLRTCC_mc,'g-*','linewidth',linewide1,'markersize',mkft);

% load Pd_CLGLRT4_2Kmu1lambda3s0.3o1_p.mat
% plot(SNRout,Pd_CLGLRT_mc,'k-v','linewidth',linewide1,'markersize',mkft);
% plot(SNRout,Pd_KGLRTCC_mc,'g-v','linewidth',linewide1,'markersize',mkft);
% load Pd_CLGLRT4_2Kmu1lambda3s0.4o1_p.mat
% plot(SNRout,Pd_CLGLRT_mc,'k-h','linewidth',linewide1,'markersize',mkft);
% plot(SNRout,Pd_KGLRTCC_mc,'g-h','linewidth',linewide1,'markersize',mkft);

load Pd_CLGLRT4_2Kmu1lambda3s0.5o1_p.mat
plot(SNRout,Pd_CLGLRT_mc,'k-<','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-<','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-<','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-<','linewidth',linewide1,'markersize',mkft);
h_leg = legend('GLC-GLRT,\sigma^2=0.1','1S-GLRT with CC,\sigma^2=0.1',...
       'GLC-GLRT,\sigma^2=0.5','1S-GLRT with CC,\sigma^2=0.5',...
       '1S-GLRT with SCM','1S-GLRT with NSCM');
% h_leg = legend('GLC-GLRT,\sigma^2=0.1','1S-GLRT with CC,\sigma^2=0.1',...
%        'GLC-GLRT,\sigma^2=0.2','1S-GLRT with CC,\sigma^2=0.2',...
%        'GLC-GLRT,\sigma^2=0.3','1S-GLRT with CC,\sigma^2=0.3',...
%        'GLC-GLRT,\sigma^2=0.4','1S-GLRT with CC,\sigma^2=0.4',...
%        'GLC-GLRT,\sigma^2=0.5','1S-GLRT with CC,\sigma^2=0.5',...
%        '1S-GLRT with SCM','1S-GLRT with NSCM');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%非高斯环境下数据量的关系%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
load Pd_CLGLRT4_1Kmu1lambda3s0.7o1_p.mat
% Pd_CLGLRT_mc(1:4)=Pd_CLGLRT_mc(1:4)*0.8;
% Pd_CLGLRT_mc(5:end)=Pd_CLGLRT_mc(1:end-4);
plot(SNRout,Pd_CLGLRT_mc,'k-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-p','linewidth',linewide1,'markersize',mkft);
load Pd_CLGLRT4_2Kmu1lambda3s0.7o1_p.mat
plot(SNRout,Pd_CLGLRT_mc,'k->','linewidth',linewide1,'markersize',mkft);%Pd_CLGLRT_mc/2+Pd_KGLRTCC_mc/2
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',linewide1,'markersize',mkft);
load Pd_CLGLRT4_3Kmu1lambda3s0.7o1_p.mat
plot(SNRout,Pd_CLGLRT_mc,'k-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-o','linewidth',linewide1,'markersize',mkft);
h_leg = legend('GLC-GLRT，K=N','1S-GLRT with CC，K=N',...
               '1S-GLRT with SCM，K=N','1S-GLRT with NSCM，K=N',...
                'GLC-GLRT，K=2N','1S-GLRT with CC，K=2N',...
                '1S-GLRT with SCM，K=2N','1S-GLRT with NSCM，K=2N',...
                 'GLC-GLRT，K=3N','1S-GLRT with CC，K=3N',...
                '1S-GLRT with SCM，K=3N','1S-GLRT with NSCM，K=3N');
% h_leg = legend('GLC-GLRT,K=8','1S-GLRT with CC,K=8',...
%                '1S-GLRT with SCM,K=8','1S-GLRT with NSCM,K=8',...
%                'GLC-GLRT,K=16','1S-GLRT with CC,K=16',...
%                '1S-GLRT with SCM,K=16','1S-GLRT with NSCM,K=16',...
%                'GLC-GLRT,K=24','1S-GLRT with CC,K=24',...
%                '1S-GLRT with SCM,K=24','1S-GLRT with NSCM,K=24');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%高斯环境下%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
hold on
load Pd_iGLCGLRT_1Kmu1lambda3s0.1o1_g.mat
% SNRout = SNRout-5;
plot(SNRout,Pd_CLGLRT_mc,'k-p','linewidth',linewide1,'markersize',mkft);%Pd_CLGLRT_mc/2+Pd_KGLRTCC_mc/2
plot(SNRout,Pd_KGLRTCC_mc,'g-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-p','linewidth',linewide1,'markersize',mkft);
load Pd_iGLCGLRT_2Kmu1lambda3s0.1o1_g.mat
% SNRout = SNRout-5;
plot(SNRout,Pd_CLGLRT_mc,'k->','linewidth',linewide1,'markersize',mkft);%Pd_CLGLRT_mc/2+Pd_KGLRTCC_mc/2
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',linewide1,'markersize',mkft);
load Pd_iGLCGLRT_3Kmu1lambda3s0.1o1_g.mat
% SNRout = SNRout-5;
plot(SNRout,Pd_CLGLRT_mc,'k-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-o','linewidth',linewide1,'markersize',mkft);
h_leg = legend('GLC-GLRT，K=N','1S-GLRT with CC，K=N',...
               '1S-GLRT with SCM，K=N','1S-GLRT with NSCM，K=N',...
                'GLC-GLRT，K=2N','1S-GLRT with CC，K=2N',...
                '1S-GLRT with SCM，K=2N','1S-GLRT with NSCM，K=2N',...
                 'GLC-GLRT，K=3N','1S-GLRT with CC，K=3N',...
                '1S-GLRT with SCM，K=3N','1S-GLRT with NSCM，K=3N');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% axis([-5,15,0,1])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%不同lambda下的检测概率%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
hold on
load Pd_CLGLRT4_2_2Kmu1lambda1s0.1o1_p.mat
plot(SNRout,Pd_CLGLRT_mc,'k-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-p','linewidth',linewide1,'markersize',mkft);
% load Pd_CLGLRT2_2Kmu1lambda1.2s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc,'k-o','linewidth',linewide1,'markersize',mkft);
% plot(SNRout_real,Pd_KGLRT_mc,'g-o','linewidth',linewide1,'markersize',mkft);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-o','linewidth',linewide1,'markersize',mkft);
% plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'b-o','linewidth',linewide1,'markersize',mkft);
% load Pd_CLGLRT4_2Kmu1lambda2s0.1o1_p.mat
% plot(SNRout(1:26),Pd_CLGLRT_mc(1:26),'k->','linewidth',linewide1,'markersize',mkft);
% plot(SNRout(1:26),Pd_KGLRTCC_mc(1:26),'g->','linewidth',linewide1,'markersize',mkft);
% plot(SNRout(1:26),Pd_KGLRT_mc(1:26),'b->','linewidth',linewide1,'markersize',mkft);
% plot(SNRout(1:26),Pd_KGLRTNSCM_mc(1:26),'c->','linewidth',linewide1,'markersize',mkft);
load Pd_CLGLRT4_2_2Kmu1lambda3s0.1o1_p.mat
plot(SNRout,Pd_CLGLRT_mc,'k->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',linewide1,'markersize',mkft);
% load Pd_CLGLRT2_2Kmu1lambda4s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'k-v','linewidth',linewide1,'markersize',mkft);
% plot(SNRout_real,Pd_KGLRT_mc(1:L),'g-v','linewidth',linewide1,'markersize',mkft);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-v','linewidth',linewide1,'markersize',mkft);
% plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'b-v','linewidth',linewide1,'markersize',mkft);
load Pd_CLGLRT4_2_2Kmu1lambda5s0.1o1_p.mat
plot(SNRout,Pd_CLGLRT_mc,'k-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-o','linewidth',linewide1,'markersize',mkft);
h_leg = legend...
       ('\lambda=1,GLC-GLRT',        '\lambda=1,1S-GLRT with CC ',...
       '\lambda=1,1S-GLRT with SCM ','\lambda=1,1S-GLRT with NSCM ',...
       '\lambda=3,GLC-GLRT',        '\lambda=3,1S-GLRT with CC ',...
       '\lambda=3,1S-GLRT with SCM ','\lambda=3,1S-GLRT with NSCM ',...
       '\lambda=5,GLC-GLRT',        '\lambda=5,1S-GLRT with CC ',...
       '\lambda=5,1S-GLRT with SCM ','\lambda=5,1S-GLRT with NSCM ');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on 
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%非高斯环境下数据量的关系2%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
load Pd_iGLCGLRT_1Kmu1lambda3s0.1o1_p.mat
% Pd_CLGLRT_mc(1:4)=Pd_CLGLRT_mc(1:4)*0.8;
% Pd_CLGLRT_mc(5:end)=Pd_CLGLRT_mc(1:end-4);
plot(SNRout,Pd_CLGLRT_mc,'k-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-p','linewidth',linewide1,'markersize',mkft);
load Pd_iGLCGLRT_2Kmu1lambda3s0.1o1_p.mat
plot(SNRout,Pd_CLGLRT_mc,'k->','linewidth',linewide1,'markersize',mkft);%Pd_CLGLRT_mc/2+Pd_KGLRTCC_mc/2
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',linewide1,'markersize',mkft);
load Pd_iGLCGLRT_3Kmu1lambda3s0.1o1_p.mat
plot(SNRout,Pd_CLGLRT_mc,'k-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-o','linewidth',linewide1,'markersize',mkft);
h_leg = legend('GLC-GLRT，K=N','1S-GLRT with CC，K=N',...
               '1S-GLRT with SCM，K=N','1S-GLRT with NSCM，K=N',...
                'GLC-GLRT，K=2N','1S-GLRT with CC，K=2N',...
                '1S-GLRT with SCM，K=2N','1S-GLRT with NSCM，K=2N',...
                 'GLC-GLRT，K=3N','1S-GLRT with CC，K=3N',...
                '1S-GLRT with SCM，K=3N','1S-GLRT with NSCM，K=3N');
% h_leg = legend('GLC-GLRT,K=8','1S-GLRT with CC,K=8',...
%                '1S-GLRT with SCM,K=8','1S-GLRT with NSCM,K=8',...
%                'GLC-GLRT,K=16','1S-GLRT with CC,K=16',...
%                '1S-GLRT with SCM,K=16','1S-GLRT with NSCM,K=16',...
%                'GLC-GLRT,K=24','1S-GLRT with CC,K=24',...
%                '1S-GLRT with SCM,K=24','1S-GLRT with NSCM,K=24');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
