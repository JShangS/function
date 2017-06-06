%%对高速高机动目标慢时间用一下LVD
clc
clear 
close all
%%%%%%%%%%%%信号参数%%%%%%%%%%%%%%%%%%%%%%%%%
fc = 1e9;
B = 1e6;
Fs = 5*B;
Ts = 1/Fs;
C = 3e8;
delt_R = C / 2 / Fs;
Tao = 128e-6;
mu = B/Tao;
t = (-Tao/2:Ts:Tao/2-Ts);
L_t =length(t);
lamda = C/(fc);
PRF = 301;
Beishu = 1;
Tr = 1/PRF;
pulse_M = Beishu*PRF;
%%%%%%%%%%%%%%%%%%%目标参数%%%%%%%%%%%%%%%%%%%%%%%%%%%
R0 = 50*delt_R;
mu_Lvd = -50;
f0_Lvd =10000;
V0 = f0_Lvd*lamda/2;
a_T = mu_Lvd*lamda/2;
jerk = 0;
%%%%%%%%%%%%%%%%%%%%目标回波和脉冲压缩%%%%%%%%%%%%%%%%%%%%%%%%
ht = exp(1j * 2 * pi * (0.5 * mu * t.^2));
ht_fft = fft(ht);
ht_fft = repmat(ht_fft,pulse_M,1);
tm = (-pulse_M/2:pulse_M/2-1)*Tr;
Rtm = R0 + V0.*tm + 0.5 * a_T * tm.^2 + 1/6 * jerk * tm.^3;
dt = repmat((2*Rtm/C)',1,L_t);
t_ML = repmat(t,pulse_M,1);
echo = exp(-1j * 2 * pi * (fc * dt + 0.5 * mu * (t-dt).^2));
echo = awgn(echo,-22);
echo_fft = fft(echo,[],2);
pc = ifft(echo_fft.*ht_fft,[],2);
pc = pc.'; %行是慢时间，列是快时间
figure()
mesh(abs(pc))
%%%%%%%%%%%%%%%%%%%%%%%%%%  MTD  %%%%%%%%%%%%%%%%%%%%%%
MTD = fft(pc,[],2);
figure()
mesh(abs(MTD))
title('MTD')
%%%%%%%%%%%
[~,index_lie] = max(abs(pc));
target_line = pc((0:pulse_M-1)*L_t + index_lie);
target_line = target_line.';
% target_line = target_line./abs(target_line); %%归一化
%%%%%%%%%%%%%%%%%%%%%%%%%%  ACCF-LVD  %%%%%%%%%%%%%
ACCF = ifftshift(ifft(fft(pc(:,1:pulse_M-1)).* conj(fft(pc(:,2:pulse_M)))),1);
figure()
mesh(abs(ACCF))
title('ACCF结果')
[~,index_lie] = max(abs(pc));
target_line_accf = ACCF(129,:);
target_line_accf = target_line_accf.';
a = 1;
q = a/Tr;
h = 1;
s1 = JS_RXC3(target_line_accf,q);
s1_fft = fftshift(fft(s1,[],1),1);
% figure()
% mesh(abs(s1_fft));
% title('参数对称瞬时自相关频谱')
% figure()
% mesh(abs(s1))
% title('参数对称瞬时自相关')

%%变尺度FT
tic
S = JS_SFT2(s1,a,h,Tr); 
toc
% figure()
% mesh(abs(S))
% title('变尺度FT')
LVD = fftshift(fft(S,[],2),2);
figure()
mesh(abs(LVD))
title('ACCF-LVD')
% figure()
% plot(abs(target_line))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RFRFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha= 0:0.001:2;
L_alpha = length(alpha);
tic
for i_a = 1:L_alpha
    i_a
    x(i_a,:) = frft(target_line,alpha(i_a));
end
toc
figure()
mesh((abs(x)));
title('FRFT结果')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% RLVD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % target_line_fft = fft(target_line);
% % figure()
% % plot(abs(target_line_fft))
% % 
% % a = 1;
% % q = a/Tr;
% % h = 1;
% % s1 = JS_RXC3(target_line,q);
% % s1_fft = fftshift(fft(s1,[],1),1);
% % figure()
% % mesh(abs(s1_fft));
% % title('参数对称瞬时自相关频谱')
% % figure()
% % mesh(abs(s1))
% % title('参数对称瞬时自相关')
% % 
% % %变尺度FT
% % tic
% % S = JS_SFT2(s1,a,h,Tr); 
% % toc
% % figure()
% % mesh(abs(S))
% % title('变尺度FT')
% % LVD = fftshift(fft(S,[],2),2);
a = 1;
q = a/Tr;
h = 1;
LVD = JS_LVD(target_line,Tr);
LVD = circshift(abs(LVD),[0,-(pulse_M-Beishu*f0_Lvd)]);
f1 = linspace(-PRF/2,PRF/2,pulse_M);
f_u = linspace(-PRF/2,PRF/2,pulse_M);
[F,Mu] = (meshgrid(f1,f_u));
figure()
mesh(-F,-Mu,abs(LVD))
% mesh(abs(LVD))
title('RLVD')
% axis([-PRF/2,PRF/2,-PRF/2,PRF/2,0,max(max(abs(LVD)))]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
