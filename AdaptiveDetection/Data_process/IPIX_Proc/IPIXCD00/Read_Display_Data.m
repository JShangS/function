
clear ; close all;clc; 

%***************************IPIX 雷达数据处理**********************************
cdfFile =  '19980205_170935_ANTSTEP.CDF';
[nc nrange nsweep ntxpol nadc cdfFileName] = ipixinfo(cdfFile,'tmp.txt');
[I, Q, meanIQ, stdIQ, inbal,adc_data]=ipixload(nc,'hh',0,'auto');
sig = I + 1j*Q;   
filename = strcat(cdfFile,'.mat');
save(filename,'sig');
% sig = sig/max(abs(sig));
load 19980205_170935_ANTSTEP.CDF.mat
sig_test=sig(:,1:end).';
figure; mesh(abs(((sig_test))));   colormap jet;

% sig = sum(sig,2);
% PRF = 1000;
% Npulse = nsweep;
% Nfft = 1024;
% tfr = tfrstft(sig,1:Npulse,Nfft,hamming(65));
% tfr = fftshift(tfr,1);
% faxis = (-Nfft/2+1:Nfft/2)/Nfft*PRF;
% taxis = (0:Npulse-1)/PRF;
% figure;imagesc(taxis,faxis,abs(tfr));  colormap jet;