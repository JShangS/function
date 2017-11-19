clc
clear 
close all
load Pd_CLGLRT_2K_mu1_lamda1_sigma0.1.mat
SNRout_real=0:1:20; % Êä³öSNR
L = length(SNRout_real)
plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-s','linewidth',2,'markersize',10);
hold on
load Pd_CLGLRT_2K_mu1_lamda1.2_sigma0.1.mat
plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-s','linewidth',2,'markersize',10);

load Pd_CLGLRT_2K_mu1_lamda2_sigma0.1.mat
plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-o','linewidth',2,'markersize',10);

load Pd_CLGLRT_2K_mu1_lamda3_sigma0.1.mat
plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-*','linewidth',2,'markersize',10);

load Pd_CLGLRT_2K_mu1_lamda4_sigma0.1.mat
plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-v','linewidth',2,'markersize',10);
load Pd_CLGLRT_2K_mu1_lamda5_sigma0.1.mat
plot(SNRout_real,Pd_CLGLRT_mc(1:L),'b-h','linewidth',2,'markersize',10);
legend('\lambda=1','\lambda=1.2','\lambda=2','\lambda=3','\lambda=4','\lambda=5')
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
grid on