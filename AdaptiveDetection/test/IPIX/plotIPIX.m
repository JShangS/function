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
load Pd_MAM_IPIX_2K_19980223_170435_.mat
% Pd_CLGLRT_mc(1:end-1)=Pd_CLGLRT_mc(2:end);
plot(SNRout,Pd_SCM_mc,'k-*','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_NSCM_mc,'g-*','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_MAM_mc,'b-*','linewidth',linewide1,'markersize',mkft);

load Pd_CLGLRT_19980223_170435_IPI_2Ks0.mat
% Pd_CLGLRT_mc(1:end-1)=Pd_CLGLRT_mc(2:end);
plot(SNRout,Pd_CLGLRT_mc,'k->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',linewide1,'markersize',mkft);

load Pd_CLGLRT_19980223_170435_IPI_3Ks0.mat
% Pd_CLGLRT_mc(1:end-1)=Pd_CLGLRT_mc(2:end);
plot(SNRout,Pd_CLGLRT_mc,'k-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-o','linewidth',linewide1,'markersize',mkft);

xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
% axis([-5,15,0,1])
grid on
box on
h_leg = legend('GLC-GLRT，K=N','1S-GLRT with CC，K=N',...
               '1S-GLRT with SCM，K=N','1S-GLRT with NSCM，K=N',...
                'GLC-GLRT，K=2N','1S-GLRT with CC，K=2N',...
                '1S-GLRT with SCM，K=2N','1S-GLRT with NSCM，K=2N',...
                 'GLC-GLRT，K=3N','1S-GLRT with CC，K=3N',...
                '1S-GLRT with SCM，K=3N','1S-GLRT with NSCM，K=3N');
set(h_leg,'Location','SouthEast')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
load Pd_CLGLRT_2_19980223_170435_IPI_3Ks0.mat
Pd_CLGLRT_mc(1:end-1)=Pd_CLGLRT_mc(2:end);
plot(SNRout,Pd_CLGLRT_mc,'k-*','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-*','linewidth',linewide1,'markersize',mkft);
load Pd_CLGLRT_2_19980223_170435_IPI_1Ks0.mat
Pd_CLGLRT_mc(1:end-1)=Pd_CLGLRT_mc(2:end);
plot(SNRout,Pd_KGLRT_mc,'b-*','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-*','linewidth',linewide1,'markersize',mkft);

load Pd_CLGLRT_2_19980223_170435_IPI_2Ks0.mat
Pd_CLGLRT_mc(1:end-1)=Pd_CLGLRT_mc(2:end);
plot(SNRout,Pd_CLGLRT_mc,'k->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',linewide1,'markersize',mkft);

load Pd_CLGLRT_2_19980223_170435_IPI_1Ks0.mat
Pd_CLGLRT_mc(1:end-1)=Pd_CLGLRT_mc(2:end);
plot(SNRout,Pd_CLGLRT_mc,'k-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-o','linewidth',linewide1,'markersize',mkft);
load Pd_CLGLRT_2_19980223_170435_IPI_3Ks0.mat
Pd_CLGLRT_mc(1:end-1)=Pd_CLGLRT_mc(2:end);
plot(SNRout,Pd_KGLRT_mc,'b-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-o','linewidth',linewide1,'markersize',mkft);

xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
% axis([-5,15,0,1])
grid on
box on
h_leg = legend('GLC-GLRT，K=N','1S-GLRT with CC，K=N',...
               '1S-GLRT with SCM，K=N','1S-GLRT with NSCM，K=N',...
                'GLC-GLRT，K=2N','1S-GLRT with CC，K=2N',...
                '1S-GLRT with SCM，K=2N','1S-GLRT with NSCM，K=2N',...
                 'GLC-GLRT，K=3N','1S-GLRT with CC，K=3N',...
                '1S-GLRT with SCM，K=3N','1S-GLRT with NSCM，K=3N');
set(h_leg,'Location','SouthEast')