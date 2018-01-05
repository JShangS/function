%%计算LVD

clc
clear
close all
%信号%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc1 = 100;%载频
mu1 = 50;%调频斜率
fc2 = 50;%载频
mu2 = -40;%调频斜率
fc3 = -70;%载频
mu3 = 250;%调频斜率
B = 1e6;%带宽
Fs = 512;%采样频率
Ts = 1/Fs;
t = (-Fs/2:Fs/2-1)*Ts;
tao = (-Fs/2:Fs/2-1)*Ts;
L = length(t);%%快时间频点数
L_tao = length(tao);
f1 = linspace(-Fs/4,Fs/4-1,L);
f_u = linspace(-Fs/2,Fs/2-1,L_tao);
A1 = 1;
A2 = 1;
A3 = 1;
s1 = A1 * exp(1j*2*pi*(fc1 * t + 0.5*mu1*t.^2))';
% angle_s1 = angle(s1);
s2 = A2 * exp(1j*2*pi*(fc2 * t + 0.5*mu2*t.^2))';
s3 = A3 * exp(1j*2*pi*(fc3 * t + 0.5*mu3*t.^2))';
s = s1+s2+s3;
s = awgn(s,10);
% s = hilbert(s);
figure()
plot(real(s))
title('信号')
s_fft = (fft(s));
figure()
plot(f1,abs(s_fft))
title('信号频谱')
%%%尺度变换参数
a = 1;
q = a/Ts;
h = 1;
s1 = JS_RXC3(s,q); %%变尺度瞬时自相关
%%%s1=  [(tao1,t1),(tao1,t2),(tao1,t3)
%%%      (tao2,t1),(tao2,t2),(tao2,t3)
%%%      (tao3,t1),(tao3,t2),(tao3,t3)]
s1_fft = fftshift(fft(s1,[],1),1);
figure()
mesh(abs(s1_fft));
% hold off
title('参数对称瞬时自相关频谱')
figure()
mesh(abs(s1))
title('参数对称瞬时自相关')
%%变尺度FT
S = JS_SFT2(s1,a,h,Ts); 
figure()
mesh(abs(S))
title('变尺度FT')
LVD = fftshift(fft(S,[],2),2);
[F,Mu] = (meshgrid(f1,f_u));
figure()
mesh(-Mu,-F,abs(LVD))
title('LV分布结果')
ylabel('中心频率(Hz)')
xlabel('调频率(Hz/s)')
axis([min(f_u),max(f_u),min(f1),max(f1),0,max(max(abs(LVD)))])

 