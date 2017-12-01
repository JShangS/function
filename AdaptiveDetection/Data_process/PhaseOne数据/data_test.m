clc 
clear
close all
load H067038_3iqHH_H067037_2iqVV.mat
N = 16;
Range = 30;
pulse = 1*N:2*N-1;
pulse_num = length(pulse);
Zhh_test = Zhh(pulse,:);
figure(1)
imagesc(abs(Zhh))
[X,Y]=meshgrid(1:76,linspace(-1,1,pulse_num));
MTD = abs(fftshift(fft(Zhh_test,[],1)));
figure(2)
mesh(X,Y,MTD);
R_KA = zeros(N,N);
L_KA = 1000;
for i = 1:L_KA
    X_KA = Zhh((i-1)*16+1:i*16,Range);
    R_KA = R_KA+X_KA*X_KA'/76/L_KA;
end
% X_KA = Zhh(40:end,Rang);
% L_KA = length(X_KA);
figure(3)
mesh(abs(R_KA))
figure()
R_test = Zhh_test(:,Range)*Zhh_test(:,Range)'/76;
mesh(abs(R_test))
% R_KA = X_KA*X_KA'/L_KA;

