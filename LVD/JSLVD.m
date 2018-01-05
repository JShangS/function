%%����LVD

clc
clear
close all
%�ź�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc1 = 100;%��Ƶ
mu1 = 50;%��Ƶб��
fc2 = 50;%��Ƶ
mu2 = -40;%��Ƶб��
fc3 = -70;%��Ƶ
mu3 = 250;%��Ƶб��
B = 1e6;%����
Fs = 512;%����Ƶ��
Ts = 1/Fs;
t = (-Fs/2:Fs/2-1)*Ts;
tao = (-Fs/2:Fs/2-1)*Ts;
L = length(t);%%��ʱ��Ƶ����
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
title('�ź�')
s_fft = (fft(s));
figure()
plot(f1,abs(s_fft))
title('�ź�Ƶ��')
%%%�߶ȱ任����
a = 1;
q = a/Ts;
h = 1;
s1 = JS_RXC3(s,q); %%��߶�˲ʱ�����
%%%s1=  [(tao1,t1),(tao1,t2),(tao1,t3)
%%%      (tao2,t1),(tao2,t2),(tao2,t3)
%%%      (tao3,t1),(tao3,t2),(tao3,t3)]
s1_fft = fftshift(fft(s1,[],1),1);
figure()
mesh(abs(s1_fft));
% hold off
title('�����Գ�˲ʱ�����Ƶ��')
figure()
mesh(abs(s1))
title('�����Գ�˲ʱ�����')
%%��߶�FT
S = JS_SFT2(s1,a,h,Ts); 
figure()
mesh(abs(S))
title('��߶�FT')
LVD = fftshift(fft(S,[],2),2);
[F,Mu] = (meshgrid(f1,f_u));
figure()
mesh(-Mu,-F,abs(LVD))
title('LV�ֲ����')
ylabel('����Ƶ��(Hz)')
xlabel('��Ƶ��(Hz/s)')
axis([min(f_u),max(f_u),min(f1),max(f1),0,max(max(abs(LVD)))])

 