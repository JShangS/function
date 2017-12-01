clc
clear 
close all
load 19980205_170935_ANTSTEP.mat
N = 8;
Range = 14;
[M,L] = size(sig_test);
figure(1)
imagesc(abs(sig_test))
sig_test_t = sig_test(2:N+1,:);
[X,Y]=meshgrid(1:L,linspace(-0.5,0.5,M));
MTD = abs(fftshift(fft(sig_test,[],1)));
figure(2)
mesh(X,Y,MTD);


R_KA = zeros(N,N);
L_KA = 7500;
for i = 1:L_KA
    X_KA = sig_test((i-1)*N+1:i*N,Range);
    R_KA = R_KA+X_KA*X_KA'/L/L_KA;
end
figure(3)
mesh(abs(R_KA))
save('19980205_170935_IPIX','R_KA','Range','N','sig_test');
