clear all
close all
clc


ncid=netcdf.open('19931107_135603_starea.cdf');
[ I , Q ]=ipixload(ncid , 'vv' , 2 , 'auto');
subplot(2,3,1)
plot(I)
xlabel('Time [ms]')
ylabel('Sea clutter returns')
axis tight

[f, x] = ksdensity(I,-10:0.01:10);
subplot(2,3,2)
plot(x, f)
xlabel('Sea clutter returns')
ylabel('Density')
ylim([0 1])
[f, x] = ksdensity(normrnd(0,var(I),1,1e5),-10:0.01:10);
subplot(2,3,3)
plot(x, f)
xlabel('Gaussian')
ylabel('Density')
ylim([0 1])