% ===========================================================================================%
%            �ó������8�������źŵ���ѹ��MTI/MTD                                           %
% ===========================================================================================%
close all; %�ر�����ͼ��
clear  %������б���
clc;
% ===================================================================================%
%                                    �״����                                       %
% ===================================================================================%
C=3.0e8;  %����(m/s)
RF=3.140e9;  %�״���Ƶ
Lambda=C/RF;%�״﹤������
PulseNumber=8;   %�ز������� 
BandWidth=2.0e6;  %�����źŴ���
TimeWidth=42.0e-6; %�����ź�ʱ��
PRT=120e-6;   % �״﷢�������ظ�����(s),120us��Ӧ1/2*120*300=18000��
PRF=1/PRT;
Fs=2.0e6;  %����Ƶ��
NoisePower=-12;%(dB);%�������ʣ�Ŀ��Ϊ0dB��
% ---------------------------------------------------------------%
SampleNumber=fix(Fs*PRT);%����һ���������ڵĲ���������
TotalNumber=SampleNumber*PulseNumber;%�ܵĲ���������
BlindNumber=fix(Fs*TimeWidth);%����һ���������ڵ�ä��-�ڵ���������
%===================================================================================%
%                                    Ŀ�����                                       %
%===================================================================================%
 TargetNumber=3;%Ŀ�����
SigPower(1:TargetNumber)=[1 1 1];%Ŀ�깦��,������
 TargetDistance(1:TargetNumber)=[3000 8025 10125];%Ŀ�����,��λm
 DelayNumber(1:TargetNumber)=fix(Fs*2*TargetDistance(1:TargetNumber)/C);% ��Ŀ����뻻��ɲ����㣨�����ţ�
TargetVelocity(1:TargetNumber)=[0 50 450];%Ŀ�꾶���ٶ� ��λm/s
TargetFd(1:TargetNumber)=2*TargetVelocity(1:TargetNumber)/Lambda; %����Ŀ��ಷ��
%====================================================================================%
%                                   �������Ե�Ƶ�ź�                                     %
%====================================================================================%
 number=fix(Fs*TimeWidth);%�ز��Ĳ�������=��ѹϵ������=��̬����Ŀ+1
if rem(number,2)~=0
   number=number+1;
end   
for i=-fix(number/2):fix(number/2)-1
  Chirp(i+fix(number/2)+1)=exp(j*(pi*(BandWidth/TimeWidth)*(i/Fs)^2));
  
end
coeff=conj(fliplr(Chirp));%Flip matrix left to right; returns the complex
                          %conjugate of the elements of Z

 
%-------------------------����Ŀ��ز���------------------------%
SignalAll=zeros(1,TotalNumber);%����������ź�,����0
for k=1:TargetNumber% ���β�������Ŀ��
   SignalTemp=zeros(1,SampleNumber);% һ������
   SignalTemp(DelayNumber(k)+1:DelayNumber(k)+number)=sqrt(SigPower(k))*Chirp;%һ�������1��Ŀ�꣨δ�Ӷ������ٶȣ�
   
   Signal=zeros(1,TotalNumber);
   for i=1:PulseNumber
      Signal((i-1)*SampleNumber+1:i*SampleNumber)=SignalTemp;
   end
   FreqMove=exp(j*2*pi*TargetFd(k)*(0:TotalNumber-1)/Fs);%Ŀ��Ķ������ٶ�*ʱ��=Ŀ��Ķ���������
   Signal=Signal.*FreqMove;
   SignalAll=SignalAll+Signal;
end
figure(2);
subplot(2,1,1);plot(real(SignalAll),'r-');title('Ŀ���źŵ�ʵ��');
subplot(2,1,2);plot(imag(SignalAll));title('Ŀ���źŵ��鲿');
grid on;zoom on;

%====================================================================================%
%                                   ����ϵͳ�����ź�                                  %
%====================================================================================%
SystemNoise=normrnd(0,10^(NoisePower/10),1,TotalNumber)+j*normrnd(0,10^(NoisePower/10),1,TotalNumber);
%====================================================================================%
%                                   �ܵĻز��ź�                                     %
%====================================================================================%
Echo=SignalAll+SystemNoise;% +SeaClutter+TerraClutter;
for i=1:PulseNumber   %�ڽ��ջ�������,���յĻز�Ϊ0
      Echo((i-1)*SampleNumber+1:(i-1)*SampleNumber+number)=0;
end
figure(3);
subplot(2,1,1);plot(real(Echo));title('�ܻز��źŵ�ʵ��,������Ϊ0');
subplot(2,1,2);plot(imag(Echo));title('�ܻز��źŵ��鲿,������Ϊ0');
%================================ʱ����ѹ=================================%
pc_time0=conv(Echo,coeff);
figure(4);plot(abs(pc_time0));title('ʱ����ѹ����ķ���,����̬��');
pc_time1=pc_time0(number:TotalNumber+number-1);%ȥ����̬�� number-1��
% ================================Ƶ����ѹ=================================%
Echo_fft=fft(Echo,2048);%��Ӧ����TotalNumber+number-1��FFT,��Ϊ����������ٶ�,������2048���FFT
coeff_fft=fft(coeff,2048);
pc_fft=Echo_fft.*coeff_fft;
pc_freq0=ifft(pc_fft);
figure(5);subplot(2,1,1);plot(abs(pc_freq0));title('Ƶ����ѹ����ķ���,��ǰ��̬��ͺ���̬��');
hold on;
subplot(2,1,2);plot(abs(pc_time0(1:TotalNumber+number-1)-pc_freq0(1:TotalNumber+number-1)),'r');title('��ɫΪʱ��Ƶ����ѹ�Ĳ��');
pc_freq1=pc_freq0(number:TotalNumber+number-1);%ȥ����̬�� number-1��,����������(2048-number+1-TotalNumber=45��)
% ================��������š������ź���������=================================%
for i=1:PulseNumber
      pc(i,1:SampleNumber)=pc_freq1((i-1)*SampleNumber+1:i*SampleNumber);
end
figure(6);
plot(abs(pc(1,:)));
% mesh(abs(pc));title('��ѹ ���');
% ================MTI����Ŀ����ʾ��,������ֹĿ��͵���Ŀ��---�������Ӳ�=================================%
mti = pc(2:end,:)-pc(1:end-1,:);
% for i=1:PulseNumber-1  %��������������һ������
%    mti(i,:)=pc(i+1,:)-pc(i,:);
% end
figure(7);
mesh(abs(mti));title('MTI ���');
% plot(abs(mti(1,:)));
% ================MTD����Ŀ���⣩,���ֲ�ͬ�ٶȵ�Ŀ�꣬�в�������=================================%
mtd = fft(pc,[],1);
% mtd=zeros(PulseNumber,SampleNumber);
% for i=1:SampleNumber
%    buff(1:PulseNumber)=pc(1:PulseNumber,i);
%    buff_fft=fft(buff);
%    mtd(1:PulseNumber,i)=buff_fft(1:PulseNumber)';
% end
%  figure(7);plot(abs(mtd));title('MTD  ���');
  figure(8);mesh(abs(mtd));title('MTD  ���');