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
%%%%%%%%%%%%%%%%%%%%%LV参数设置%%%%%%%%%%%%%%%%%%%%%%%%%
a = 1;
q = a/Tr;
h = 1;
%%%%%%%%%%%%%%%%%%%目标参数%%%%%%%%%%%%%%%%%%%%%%%%%%%
R01 = 300*delt_R;  %相对初始距离
R02 = 320*delt_R;
mu_Lvd1 = 50;  %调频率
f0_Lvd1 =100;   %多普勒
mu_Lvd2 = -50;
f0_Lvd2 =120;
V01 = f0_Lvd1 * lamda/2;  %多普勒对应的速度
a_T1 = mu_Lvd1 * lamda/2; %调频率对应的加速度
V02 = f0_Lvd2 * lamda/2;
a_T2 = mu_Lvd2 * lamda/2;
A1 = 1;
A2 = 1;
%%%%%%%%%%%%%%%%%%%%目标回波和脉冲压缩%%%%%%%%%%%%%%%%%%%%%%%%
ht = exp(1j * 2 * pi * (0.5 * mu * t.^2)); 
ht_fft = fft(ht);   %脉冲压缩系数
ht_fft = repmat(ht_fft,pulse_M,1);
tm = (-pulse_M/2:pulse_M/2-1)*Tr;   %慢时间
Rtm1 = R01 + V01.*tm + 0.5 * a_T1 * tm.^2;    %目标1对应的慢时间距离走动
Rtm2 = R02 + V02.*tm + 0.5 * a_T2 * tm.^2;    %目标2对应的慢时间距离走动
dt1 = repmat((2*Rtm1/C)',1,L_t);
dt2 = repmat((2*Rtm2/C)',1,L_t);    %目标2延迟
t_ML = repmat(t,pulse_M,1);
SNR = -10;  %信噪比
%%%%%%%%%%%%%%%%%%%%%%%%%目标1回波
echo1 = A1 * exp(-1j * 2 * pi * (fc * dt1 + 0.5 * mu * (t-dt1).^2));
echo_fft1 = fft(echo1,[],2);
%%%%%%%%%%%%%%%%%%%%%%%%%目标2回波
echo2 = A2 * exp(-1j * 2 * pi * (fc * dt2 + 0.5 * mu * (t-dt2).^2));
echo_fft2 = fft(echo2,[],2);
% %%%%%%%%%%%%%%%%%%%%%%%%%目标1回波脉压
% pc1 = ifft(echo_fft1.*ht_fft,[],2);
% pc1 = pc1.'; %行是慢时间，列是快时间
% %%%%%%%%%%%%%%%%%%%%%%%%目标2回波脉压
% pc2 = ifft(echo_fft2.*ht_fft,[],2);
% pc2 = pc2.'; %行是慢时间，列是快时间
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%直接找目标12的轨迹，先不用Radon%%%%%%%%%
% %%%%%%%%%%%找到目标1脉压对应的位置，先不用Radon找
% [~,index_lie1] = max(abs(pc1));
% target_line1 = pc1((0:pulse_M-1)*L_t + index_lie1);
% target_line1 = target_line1.';
% %%%%%%%%%%%找到目标2脉压对应的位置，先不用Radon找
% [~,index_lie2] = max(abs(pc2));
% target_line2 = pc2((0:pulse_M-1)*L_t + index_lie2);
% target_line2 = target_line2.';
% %%%%%%%% 两个目标轨迹叠加
% target_line = target_line1 + target_line2;
%%%%%%%%%%%%%%%%%%%%%%%%%%目标12回波叠加,脉压%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo = echo1 + echo2;
echo_fft = fft(echo,[],2);
echo = awgn(echo , SNR);
pc = ifft(echo_fft.*ht_fft,[],2);
pc = pc.'; %行是慢时间，列是快时间
figure()
mesh(abs(pc))
title('脉压结果')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Radon变换找轨迹%%%%%%%%%%%%%%%%%
target_line = zeros(pulse_M,1);
d_theta = atan(1/(pulse_M-1))/2;  %搜索步长,相参积累时间内移动了一个采样单元为步长
theta = -100*d_theta:d_theta:100*d_theta;
L_theta = length(theta);
for i = 1:L_theta
    i
    offset = round((0:pulse_M-1)*tan(theta(i))); %偏移量
    for R0 = 290:330
        offset = offset + R0;
        if min(offset) > 0 && max(offset) < L_t
            target_line = target_line + pc((0:pulse_M-1)*L_t + offset).';
        end
    end  
end
%%%%%%%%%%%%%%%%%%%  MTD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MTD = fft(pc,[],2);
figure()
mesh(abs(MTD))
title('MTD')
%%%%%%%%%%%%%%%%%%%% ACCF-LVD %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%计算ACCF
% ACCF = ifftshift(ifft(fft(pc(:,1:pulse_M-1)).* conj(fft(pc(:,2:pulse_M)))),1);
% figure()
% mesh(abs(ACCF))
% title('ACCF结果')
% %%%%计算LVD
% target_line_accf = ACCF(L_t/2+1,:);
% target_line_accf = target_line_accf.';
% ACCF_LVD = JS_LVD(target_line_accf,Tr);
% figure()
% mesh(abs(ACCF_LVD))
% title('ACCF-LVD结果')
%%%%%%%%%%%%%%%%%%%%%%% radon-LVD %%%%%%%%%%%%%%%%%%%%%%%
% Radon_LVD1 = JS_LVD(target_line1,Tr);
% Radon_LVD1 = circshift(abs(Radon_LVD1),[0,-(pulse_M-Beishu*f0_Lvd1)]);
% Radon_LVD2 = JS_LVD(target_line2,Tr);
% Radon_LVD2 = circshift(abs(Radon_LVD2),[0,-(pulse_M-Beishu*f0_Lvd2)]);
% Radon_LVD = Radon_LVD1+Radon_LVD2;
Radon_LVD = JS_LVD(target_line,Tr);
Radon_LVD = circshift(abs(Radon_LVD),[0,-(pulse_M-Beishu*f0_Lvd2)]);
% Radon_LVD2 = circshift(abs(Radon_LVD),[0,(pulse_M-Beishu*f0_Lvd2)]);
f1 = linspace(-PRF/2,PRF/2,pulse_M);
f_u = linspace(-PRF/2,PRF/2,pulse_M);
[F,Mu] = (meshgrid(f1,f_u));
figure()
% mesh(abs(Radon_LVD))
mesh(F,Mu,(abs(Radon_LVD)));%fftshift
title('Radon-LVD结果')
xlabel('f0\ ')
ylabel('mu')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RFRFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 0:0.001:2;
p=p.';
S = PRF / (Tr * pulse_M); %量纲归一化
dx = sqrt((Tr * pulse_M) * PRF);
u = ((-pulse_M/2):(pulse_M/2-1)) / (dx);
L_alpha = length(p);
tic
for i_a = 1:L_alpha
    i_a
    x(i_a,:) = frft(target_line,p(i_a));
end
toc
figure()
[X,Y] = meshgrid(u , p);
mesh(X,Y,(abs(x)));
title('FRFT结果')