%*****Time Doppler plot**********
clear; close all; clc
cdfFile = '#17.cdf';   rangebin = 7;  xpix=250;  ypix=80;

[nc nrange nsweep ntxpol nadc cdfFileName] = ipixinfo(cdfFile,'tmp.txt');
[I, Q, meanIQ, stdIQ, inbal,adc_data]=ipixload(nc,'vv',0,'auto');
for i=1:size(I,1)
    I(i,:)=smooth(I(i,:),3);
    Q(i,:)=smooth(Q(i,:),3);
end

y = I + 1j*Q;

figure; imagesc(abs(y'));  xlabel('Range Cell');  ylabel('Number of pulses');

[logTD,time,doppl]=tdoppl(nc,I,Q,rangebin);
image(time,doppl,logTD);  xlabel('Times');  ylabel('Doppler');

% trslice(nc,I,Q,300,100);

netcdf.close(nc)