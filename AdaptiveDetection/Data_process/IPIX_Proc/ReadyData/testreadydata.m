close all
clear all
clc


% hi=textread('lo.dat');

ncid=netcdf.open('F:\实测数据\IPIX Data（范）\IPIX DATA\starea\19931107_135603_starea.cdf');
[ I , Q , mnIQ , sdIQ] = ipixload( ncid , 'hh' , 9 , 'auto' );

[logTD,time,doppl]=tdoppl(ncid,I,Q,9);


% y = abs( I + 1i*Q ) * sqrt( prod( sdIQ ) );

y=I+1i*Q;
clear I Q;

% sig = y - mean(y);
sig = y ;

%****************matlab自带函数*************************************
s = spectrogram(sig,200,100,256);
 PRF = 1000;  Npulse = length(sig);
s = fftshift(s,1);
taxis = (0:size(s,2)-1)/PRF;
faxis = (-Npulse/2+1:Npulse/2)/Npulse*PRF; %
figure;imagesc(taxis,faxis,abs(s));
% figure;imagesc(taxis,faxis,circshift(abs(s),-size(s,1)/2));

% %****************短时傅里叶变换*************************************
% PRF = 1000;   h1 = window('hanning',55);
% [tfr,t,f] = tfrstft(y',300,200,h1);
% Npulse = length(f);
% taxis = (0:size(f)-1)/PRF;
% faxis = (-Npulse/2+1:Npulse/2)/Npulse*PRF; 
% figure;imagesc(taxis,faxis,circshift(abs(tfr),-size(tfr,1)/2));