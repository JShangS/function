clear all
close all
clc

ncid=netcdf.open('19931107_135603_starea.cdf');
length=1e2;

[ I_hh , Q_hh ]=ipixload(ncid , 'hh' , 6 , 'auto');
[ I_hv , Q_hv ]=ipixload(ncid , 'hv' , 6 , 'auto');
[ I_vv , Q_vv ]=ipixload(ncid , 'vv' , 6 , 'auto');
[ I_vh , Q_vh ]=ipixload(ncid , 'vh' , 6 , 'auto');

acf = autocorr(I_hh,length)';
% acf = acf(1:end-1);
acf=abs(acf);
temp=fliplr(acf);
temp=temp(1:end-1);
a2=[temp acf];
subplot(2,2,1)
plot(-length:1:length,a2)
ylim([0 1])
xlim([-length length])
xlabel('Time [ms]')
title('HH channel')

acf = autocorr(I_hv,length)';
acf=abs(acf);
temp=fliplr(acf);
temp=temp(1:end-1);
a2=[temp acf];
subplot(2,2,2)
plot(-length:1:length,a2)
ylim([0 1])
xlim([-length length])
xlabel('Time [ms]')
title('HV channel')

acf = autocorr(I_vv,length)';
acf=abs(acf);
temp=fliplr(acf);
temp=temp(1:end-1);
a2=[temp acf];
subplot(2,2,3)
plot(-length:1:length,a2)
ylim([0 1])
xlim([-length length])
xlabel('Time [ms]')
title('VV channel')

acf = autocorr(I_vh,length)';
acf=abs(acf);
temp=fliplr(acf);
temp=temp(1:end-1);
a2=[temp acf];
subplot(2,2,4)
plot(-length:1:length,a2)
ylim([0 1])
xlim([-length length])
xlabel('Time [ms]')
title('HV channel')


