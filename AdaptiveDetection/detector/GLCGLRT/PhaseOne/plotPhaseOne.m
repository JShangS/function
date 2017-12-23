clc
clear 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%实测数据IPIX%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labeltsize=20;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=2;
mkft = 10;

figure
hold on
load Pd_CLGLRT4_PhaseOne_1Ks0.mat
plot(SNRout,Pd_CLGLRT_mc,'k-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-o','linewidth',linewide1,'markersize',mkft);

% load Pd_CLGLRT4_PhaseOne_1.25Ks0_g.mat
% plot(SNRout,Pd_CLGLRT_mc,'k-s','linewidth',linewide1,'markersize',10);
% plot(SNRout,Pd_KGLRTCC_mc,'g-s','linewidth',linewide1,'markersize',10);
% plot(SNRout,Pd_KGLRT_mc,'b-s','linewidth',linewide1,'markersize',10);
% plot(SNRout,Pd_KGLRTNSCM_mc,'c-s','linewidth',linewide1,'markersize',10);
% 
% load Pd_CLGLRT4_PhaseOne_1.5Ks0_g.mat
% plot(SNRout,Pd_CLGLRT_mc,'k-p','linewidth',linewide1,'markersize',10);
% plot(SNRout,Pd_KGLRTCC_mc,'g-p','linewidth',linewide1,'markersize',10);
% plot(SNRout,Pd_KGLRT_mc,'b-p','linewidth',linewide1,'markersize',10);
% plot(SNRout,Pd_KGLRTNSCM_mc,'c-p','linewidth',linewide1,'markersize',10);

load Pd_CLGLRT4_PhaseOne_2Ks0.mat
plot(SNRout,Pd_CLGLRT_mc,'k-*','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-*','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-*','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-*','linewidth',linewide1,'markersize',mkft);

load Pd_CLGLRT4_PhaseOne_3Ks0.mat
plot(SNRout,Pd_CLGLRT_mc,'k->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',linewide1,'markersize',mkft);

xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
axis([-5,15,0,1])
grid on
box on
h_leg = legend('GLC-GLRT，K=N','1S-GLRT with CC，K=N',...
               '1S-GLRT with SCM，K=N','1S-GLRT with NSCM，K=N',...
                'GLC-GLRT，K=2N','1S-GLRT with CC，K=2N',...
                '1S-GLRT with SCM，K=2N','1S-GLRT with NSCM，K=2N',...
                 'GLC-GLRT，K=3N','1S-GLRT with CC，K=3N',...
                '1S-GLRT with SCM，K=3N','1S-GLRT with NSCM，K=3N');
set(h_leg,'Location','SouthEast')
% h_leg = legend('GLC-GLRT，K=8','1S-GLRT with CC，K=8',...
%                '1S-GLRT with SCM，K=8','1S-GLRT with NSCM，K=8',...
%                'GLC-GLRT，K=10','1S-GLRT with CC，K=10',...
%                '1S-GLRT with SCM，K=10','1S-GLRT with NSCM，K=10',...
%                'GLC-GLRT，K=12','1S-GLRT with CC，K=12',...
%                '1S-GLRT with SCM，K=12','1S-GLRT with NSCM，K=12',...
%                 'GLC-GLRT，K=16','1S-GLRT with CC，K=16',...
%                 '1S-GLRT with SCM，K=16','1S-GLRT with NSCM，K=16',...
%                 'GLC-GLRT，K=24','1S-GLRT with CC，K=24',...
%                 '1S-GLRT with SCM，K=24','1S-GLRT with NSCM，K=24');
% set(h_leg,'Location','SouthEast')