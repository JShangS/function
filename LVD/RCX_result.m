%%直接根据公式写一下变尺度瞬时自相关的结果
clc
clear
close all
fc1 = 100;%载频
fc2 = 250;
fc3 = 52;
fc4 = 80;
mu1 = 120;%调频斜率
mu2 = 50;
mu3 = 60;
mu4 = -60;
Fs = 256 ;%采样频率
Ts = 1/Fs;
t = (-Fs/2:Fs/2-1)*Ts;
L = length(t);%%快时间频点数
f = linspace(-Fs/2,Fs/2-1,L);
s1 = exp(1j*2*pi*(fc1 * t + 0.5*mu1*(t).^2))';
M = L;
m = -M/2:M/2-1;
tao = m*Ts;
a = 1;
q = a/Ts;
h = 1;
A1 = 1;
A2 = 0;
A3 = 0;
A4 = 0;
RXC_result1 = zeros(L,L);
RXC_result2 = zeros(L,L);
RXC_result3 = zeros(L,L);
RXC_result4 = zeros(L,L);
SNR = 10;  %dB
for it = 1:L% 每一列随着tao变
    RXC_result1(:,it) =  A1 * exp(1j*2*pi*(fc1*(a+tao) + mu1*t(it).*(a+2*tao)));
    RXC_result2(:,it) =  A2 * exp(1j*2*pi*(fc2*(a+tao) + mu2*t(it).*(a+2*tao)));
    RXC_result3(:,it) =  A3 * exp(1j*2*pi*(fc3*(a+tao) + mu3*t(it).*(a+2*tao)));
    RXC_result4(:,it) =  A4 * exp(1j*2*pi*(fc4*(a+tao) + mu4*t(it).*(a+2*tao)));
end
RXC_result = RXC_result1 + RXC_result2 + RXC_result3 + RXC_result4;
RXC_result = awgn(RXC_result,SNR);
% figure()
% mesh(abs(RXC_result))
RXC_result_fft = fft(RXC_result,[],1);
figure()
mesh(abs(RXC_result_fft))
%变尺度的FT
S = JS_SFT2(RXC_result,a,h,Ts);
figure()
mesh(abs(S))
%LVD
LVD = fftshift(fft(S,[],2),2);
[Mu ,F] = meshgrid(f,f);
figure()
mesh(Mu,F,abs(LVD))
axis([-Fs/2,Fs/2,-Fs/2,Fs/2,0,8e4])
% S_fft = fft(S_fft,[],2);
% figure()
% mesh(abs(S_fft))