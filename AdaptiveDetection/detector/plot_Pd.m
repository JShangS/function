clc
clear 
close all

% figure(1)
% hold on
% load Pd_CLGLRT_2Kmu1lambda3s0o1_p.mat
% SNRout_real=0:1:20; % 输出SNR
% L = length(SNRout_real);
% 
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-s','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-s','linewidth',2,'markersize',10);
% 
% load Pd_CLGLRT_2Kmu1lambda3s0.01o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b->','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r->','linewidth',2,'markersize',10);
% 
% load Pd_CLGLRT_2Kmu1lambda3s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-o','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-o','linewidth',2,'markersize',10);
% 
% load Pd_CLGLRT_2Kmu1lambda3s0.2o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-*','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-*','linewidth',2,'markersize',10);
% 
% load Pd_CLGLRT_2Kmu1lambda3s0.3o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-v','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-v','linewidth',2,'markersize',10);
% load Pd_CLGLRT_2Kmu1lambda3s0.4o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-h','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-h','linewidth',2,'markersize',10);
% 
% load Pd_CLGLRT_2Kmu1lambda3s0.5o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-<','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-<','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRT_mc(1:L),'g-<','linewidth',2,'markersize',10);
% legend('CLGLRT,\sigma^2=0','GLRTCC,\sigma^2=0',...
%     'CLGLRT,\sigma^2=0.01','GLRTCC,\sigma^2=0.01',...
%     'CLGLRT,\sigma^2=0.1','GLRTCC,\sigma^2=0.1',...
%     'CLGLRT,\sigma^2=0.2','GLRTCC,\sigma^2=0.2',...
%     'CLGLRT,\sigma^2=0.3','GLRTCC,\sigma^2=0.3',...
%     'CLGLRT,\sigma^2=0.4','GLRTCC,\sigma^2=0.4',...
%     'CLGLRT,\sigma^2=0.5','GLRTCC,\sigma^2=0.5',...
%     'KGLRT')
% xlabel('SNR/dB','FontSize',20)
% ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%非高斯环境下数据量的关系%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% hold on
% load Pd_CLGLRT_1.5Kmu1lambda3s0.1o1_p.mat
% SNRout_real=0:1:20; % 输出SNR
% L = length(SNRout_real);
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-s','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-s','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRT_mc(1:L),'g-s','linewidth',2,'markersize',10);
% load Pd_CLGLRT_2Kmu1lambda3s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b->','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r->','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRT_mc(1:L),'g->','linewidth',2,'markersize',10);
% load Pd_CLGLRT_4Kmu1lambda3s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-o','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-o','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRT_mc(1:L),'g-o','linewidth',2,'markersize',10);
% legend('CLGLRT，K=24','GLRTCC，K=24','KGLRTCC，K=24',...
%        'CLGLRT，K=32','GLRTCC，K=32','KGLRTCC，K=32',...
%        'CLGLRT，K=64','GLRTCC，K=64','KGLRTCC，K=64')
% xlabel('SNR/dB','FontSize',20)
% ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%高斯环境下%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)
% hold on
% load Pd_CLGLRT_1.5Kmu1lambda3s0.1o1_g.mat
% SNRout_real=0:1:20; % 输出SNR
% L = length(SNRout_real);
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-s','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-s','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRT_mc(1:L),'g-s','linewidth',2,'markersize',10);
% load Pd_CLGLRT_2Kmu1lambda3s0.1o1_g.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b->','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r->','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRT_mc(1:L),'g->','linewidth',2,'markersize',10);
% load Pd_CLGLRT_4Kmu1lambda3s0.1o1_g.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-o','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'r-o','linewidth',2,'markersize',10);
% plot(SNRout_real,Pd_KGLRT_mc(1:L),'g-o','linewidth',2,'markersize',10);
% legend('CLGLRT,K=24','GLRTCC,K=24','KGLRTCC,K=24',...
%        'CLGLRT,K=32','GLRTCC,K=32','KGLRTCC,K=32',...
%        'CLGLRT,K=64','GLRTCC,K=64','KGLRTCC,K=64')
% xlabel('SNR/dB','FontSize',20)
% ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%不同lambda下的检测概率%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(4)
% hold on
% load Pd_CLGLRT_2Kmu1lambda1s0.1o1_p.mat
% SNRout_real=0:1:20; % 输出SNR
% L = length(SNRout_real);
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-s','linewidth',2,'markersize',10);
% load Pd_CLGLRT_2Kmu1lambda1.2s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-o','linewidth',2,'markersize',10);
% load Pd_CLGLRT_2Kmu1lambda2s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b->','linewidth',2,'markersize',10);
% load Pd_CLGLRT_2Kmu1lambda3s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-<','linewidth',2,'markersize',10);
% load Pd_CLGLRT_2Kmu1lambda4s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-v','linewidth',2,'markersize',10);
% load Pd_CLGLRT_2Kmu1lambda5s0.1o1_p.mat
% plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-p','linewidth',2,'markersize',10);
% legend('\lambda=1','\lambda=1.2','\lambda=2','\lambda=3','\lambda=4',...
%        '\lambda=5');
% xlabel('SNR/dB','FontSize',20)
% ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% grid on 