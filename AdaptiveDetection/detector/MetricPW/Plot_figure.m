clc
clear 
close all
SNRout=-5:1:25; % Êä³öSNR
figure(1)
load 1.mat
hold on
plot(SNRout,Pd_NMF_mc,'k','linewidth',1)
plot(SNRout,Pd_SCM_mc,'k-+','linewidth',1)
plot(SNRout,Pd_NSCM_mc,'k.-','linewidth',1)
plot(SNRout,Pd_KL_mc,'k-*','linewidth',1)
plot(SNRout,Pd_CC_mc,'k-s','linewidth',1)
load 2.mat
plot(SNRout,Pd_NMF_mc,'r','linewidth',1)
plot(SNRout,Pd_SCM_mc,'r-+','linewidth',1)
plot(SNRout,Pd_NSCM_mc,'r.-','linewidth',1)
plot(SNRout,Pd_KL_mc,'r-*','linewidth',1)
plot(SNRout,Pd_CC_mc,'r-s','linewidth',1)