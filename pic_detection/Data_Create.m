clc
clear
close all
%%%产生训练和测试用的数据%%%%%
fc=1000e6;%载频%2000
B=1e6;%带宽%4
Tao=127e-6;%脉宽
Fs=2*B;%采样频率
t=-Tao/2:1/Fs:Tao/2-1/Fs;%脉冲时间
L=length(t);
mu=B/Tao;%条频率
C=3e8;
R_max=C*Tao/2;
delt_R=C/(2*Fs);%%采样距离单元
R_start1=220*delt_R;%目标1初始位置%170
lamda=C/fc;
PRF=10000;
Tr=1/PRF;
Vr_start1=100;%1初始速度
pusle_num=64;%脉冲数
PV=PRF*lamda/4;
M=pusle_num;
Pfa = 1e-6;     %虚警率
r0 = -log(Pfa); %门限
echo=zeros(M,L);%回波
for i=1:pusle_num
    Vr1(i)=Vr_start1;
    fd1(i)=2*Vr1(i)/lamda;
    delt_t1(i)=2*(R_start1+Vr1(i)*Tr*(i-1))/C;%回拨延迟
end
%%脉压系数
    ht_t=exp(-1j*2*pi*(mu/2*(t).^2));
    ht=conj(fliplr(ht_t));
    ht_fft=fftshift(fft(ht));%fftshift
Createnum = 2000;
trainLabels = zeros(Createnum/2,1);
testLabels =  zeros(Createnum/2,1);
trainData = zeros(64,64,Createnum/2);
testData = zeros(64,64,Createnum/2);
SNR=-20;
h = waitbar(0,'Please wait...');
for i_Createnum = 1:Createnum
    waitbar(i_Createnum/Createnum,h,sprintf([num2str(i_Createnum/Createnum*100),'%%']));
    if i_Createnum<=Createnum/2 %%训练样本
        if rand()>0.5 %有目标
           A1 = 1; %%1.8%1.45 
           trainLabels(i_Createnum) = 1;
        else %%没目标
           A1 = 0; %%1.8%1.45 
           trainLabels(i_Createnum) = 2;
        end
    else %%测试样本
        if rand()>0.5 %有目标
           A1 = 1; %%1.8%1.45 
           testLabels(i_Createnum-Createnum/2) = 1;
        else %%没目标
           A1 = 0; %%1.8%1.45 
           testLabels(i_Createnum-Createnum/2) = 2;
        end
    end
    for i=1:pusle_num
        echo(i,:)=A1*exp(-1j*2*pi*(mu/2*(t+delt_t1(i)).^2)+-1j*2*pi*(fc)*(delt_t1(i)));
        echo(i,:)=awgn(echo(i,:),SNR);%%加噪声
        echo_fft(i,:)=(fft(echo(i,:)));%fftshift
        pc_result(i,:)=(ifft(echo_fft(i,:).*ht_fft));%ifftshift
        pc_result_fft(i,:)=(fft(pc_result(i,:)));%fftshift%快时间FFT
    end
    MTD = abs(fft(pc_result,[],1));
%     max_MTD = max(max(MTD));
%     if A1~=0
%         MTD = MTD/max(max(MTD));
%     end
%     MTD = awgn(MTD,SNR);
%     MTD = round(MTD/max_MTD*255);
    if i_Createnum<=Createnum/2 %%训练样本
       trainData(:,:,i_Createnum) = (MTD(1:end,1:64));
    else
       testData(:,:,i_Createnum-Createnum/2) = (MTD(1:end,1:64));
    end
end
close(h)
if SNR>=0
    str=['MTDData','_',num2str(abs(SNR)),'dB_V',num2str(Vr_start1),'.mat'];
else
    str=['MTDData','_fu',num2str(abs(SNR)),'dB_V',num2str(Vr_start1),'.mat'];
end
save(str,'trainData','testData','trainLabels','testLabels')


% % % 脉压后原始回波
% figure(1), 
% mesh(abs(pc_result))
% % imagesc (abs(pc_result))
% xlabel('距离单元')
% ylabel('脉冲数')
%%MTD
% dV = PV/M;
% Vrt = -PV/2+dV:dV:PV/2;


% [X,Y] = meshgrid(1:L,Vrt);
% figure(2)
% imagesc(1:L,Vrt,abs(MTD))
% MTD = abs2dB(MTD);
% mesh(X,Y,MTD)
% xlabel('距离单元')
% ylabel('速度m/s')
% zlabel('幅度')
% MTD = uint8(MTD);
% figure(3)
% imshow(MTD(1:end,1:64))