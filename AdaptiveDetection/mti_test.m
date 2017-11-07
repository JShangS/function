% ===========================================================================================%
%            该程序完成8个脉冲信号的脉压、MTI/MTD                                           %
% ===========================================================================================%
close all; %关闭所有图形
clear  %清除所有变量
clc;
% ===================================================================================%
%                                    雷达参数                                       %
% ===================================================================================%
C=3.0e8;  %光速(m/s)
RF=3.140e9;  %雷达射频
Lambda=C/RF;%雷达工作波长
PulseNumber=8;   %回波脉冲数 
BandWidth=2.0e6;  %发射信号带宽
TimeWidth=42.0e-6; %发射信号时宽
PRT=120e-6;   % 雷达发射脉冲重复周期(s),120us对应1/2*120*300=18000米
PRF=1/PRT;
Fs=2.0e6;  %采样频率
NoisePower=-12;%(dB);%噪声功率（目标为0dB）
% ---------------------------------------------------------------%
SampleNumber=fix(Fs*PRT);%计算一个脉冲周期的采样点数；
TotalNumber=SampleNumber*PulseNumber;%总的采样点数；
BlindNumber=fix(Fs*TimeWidth);%计算一个脉冲周期的盲区-遮挡样点数；
%===================================================================================%
%                                    目标参数                                       %
%===================================================================================%
 TargetNumber=3;%目标个数
SigPower(1:TargetNumber)=[1 1 1];%目标功率,无量纲
 TargetDistance(1:TargetNumber)=[3000 8025 10125];%目标距离,单位m
 DelayNumber(1:TargetNumber)=fix(Fs*2*TargetDistance(1:TargetNumber)/C);% 把目标距离换算成采样点（距离门）
TargetVelocity(1:TargetNumber)=[0 50 450];%目标径向速度 单位m/s
TargetFd(1:TargetNumber)=2*TargetVelocity(1:TargetNumber)/Lambda; %计算目标多卜勒
%====================================================================================%
%                                   产生线性调频信号                                     %
%====================================================================================%
 number=fix(Fs*TimeWidth);%回波的采样点数=脉压系数长度=暂态点数目+1
if rem(number,2)~=0
   number=number+1;
end   
for i=-fix(number/2):fix(number/2)-1
  Chirp(i+fix(number/2)+1)=exp(j*(pi*(BandWidth/TimeWidth)*(i/Fs)^2));
  
end
coeff=conj(fliplr(Chirp));%Flip matrix left to right; returns the complex
                          %conjugate of the elements of Z

 
%-------------------------产生目标回波串------------------------%
SignalAll=zeros(1,TotalNumber);%所有脉冲的信号,先填0
for k=1:TargetNumber% 依次产生各个目标
   SignalTemp=zeros(1,SampleNumber);% 一个脉冲
   SignalTemp(DelayNumber(k)+1:DelayNumber(k)+number)=sqrt(SigPower(k))*Chirp;%一个脉冲的1个目标（未加多普勒速度）
   
   Signal=zeros(1,TotalNumber);
   for i=1:PulseNumber
      Signal((i-1)*SampleNumber+1:i*SampleNumber)=SignalTemp;
   end
   FreqMove=exp(j*2*pi*TargetFd(k)*(0:TotalNumber-1)/Fs);%目标的多普勒速度*时间=目标的多普勒相移
   Signal=Signal.*FreqMove;
   SignalAll=SignalAll+Signal;
end
figure(2);
subplot(2,1,1);plot(real(SignalAll),'r-');title('目标信号的实部');
subplot(2,1,2);plot(imag(SignalAll));title('目标信号的虚部');
grid on;zoom on;

%====================================================================================%
%                                   产生系统噪声信号                                  %
%====================================================================================%
SystemNoise=normrnd(0,10^(NoisePower/10),1,TotalNumber)+j*normrnd(0,10^(NoisePower/10),1,TotalNumber);
%====================================================================================%
%                                   总的回波信号                                     %
%====================================================================================%
Echo=SignalAll+SystemNoise;% +SeaClutter+TerraClutter;
for i=1:PulseNumber   %在接收机闭锁期,接收的回波为0
      Echo((i-1)*SampleNumber+1:(i-1)*SampleNumber+number)=0;
end
figure(3);
subplot(2,1,1);plot(real(Echo));title('总回波信号的实部,闭锁期为0');
subplot(2,1,2);plot(imag(Echo));title('总回波信号的虚部,闭锁期为0');
%================================时域脉压=================================%
pc_time0=conv(Echo,coeff);
figure(4);plot(abs(pc_time0));title('时域脉压结果的幅度,有暂态点');
pc_time1=pc_time0(number:TotalNumber+number-1);%去掉暂态点 number-1个
% ================================频域脉压=================================%
Echo_fft=fft(Echo,2048);%理应进行TotalNumber+number-1点FFT,但为了提高运算速度,进行了2048点的FFT
coeff_fft=fft(coeff,2048);
pc_fft=Echo_fft.*coeff_fft;
pc_freq0=ifft(pc_fft);
figure(5);subplot(2,1,1);plot(abs(pc_freq0));title('频域脉压结果的幅度,有前暂态点和后暂态点');
hold on;
subplot(2,1,2);plot(abs(pc_time0(1:TotalNumber+number-1)-pc_freq0(1:TotalNumber+number-1)),'r');title('红色为时域频域脉压的差别');
pc_freq1=pc_freq0(number:TotalNumber+number-1);%去掉暂态点 number-1个,后填充点若干(2048-number+1-TotalNumber=45个)
% ================按照脉冲号、距离门号重排数据=================================%
for i=1:PulseNumber
      pc(i,1:SampleNumber)=pc_freq1((i-1)*SampleNumber+1:i*SampleNumber);
end
figure(6);
plot(abs(pc(1,:)));
% mesh(abs(pc));title('脉压 结果');
% ================MTI（动目标显示）,对消静止目标和低速目标---可抑制杂波=================================%
mti = pc(2:end,:)-pc(1:end-1,:);
% for i=1:PulseNumber-1  %滑动对消，少了一个脉冲
%    mti(i,:)=pc(i+1,:)-pc(i,:);
% end
figure(7);
mesh(abs(mti));title('MTI 结果');
% plot(abs(mti(1,:)));
% ================MTD（动目标检测）,区分不同速度的目标，有测速作用=================================%
mtd = fft(pc,[],1);
% mtd=zeros(PulseNumber,SampleNumber);
% for i=1:SampleNumber
%    buff(1:PulseNumber)=pc(1:PulseNumber,i);
%    buff_fft=fft(buff);
%    mtd(1:PulseNumber,i)=buff_fft(1:PulseNumber)';
% end
%  figure(7);plot(abs(mtd));title('MTD  结果');
  figure(8);mesh(abs(mtd));title('MTD  结果');