clc
clear
close all
%%%MTD和AMF的检测概率比较，MTD结果进行了CA-CFAR检测，AMF直接和门限值进行比较得到检测结果。
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
R_start1=L/2*delt_R;%目标1初始位置%170
lamda=C/fc;
PRF=10000;
Tr=1/PRF;
Vr_start1=290;%1初始速度
fd = 2*Vr_start1/lamda;
pusle_num=12;%脉冲数
PV=PRF*lamda/4;
M=pusle_num;
dV = PV/M;
Vrt = -PV/2+dV:dV:PV/2;
echo=zeros(M,L);%回波
Pfa = 1e-3;    % 虚警率
r0 = -log(Pfa); % AMF门限<A new CFAR detection test for radar>给出的
S = exp(-1j*2*pi*fd*(0:M-1).'*Tr);  %导向矢量
for i=1:pusle_num
    Vr1(i)=Vr_start1;
    fd1(i)=2*Vr1(i)/lamda;
    delt_t1(i)=2*(R_start1+Vr1(i)*Tr*(i-1))/C;%回拨延迟
end
%%脉压系数
    ht_t=exp(-1j*2*pi*(mu/2*(t).^2));
    ht=conj(fliplr(ht_t));
    ht_fft=fftshift(fft(ht));%fftshift
A1=1; 
SNR = 100;
for i=1:pusle_num
    echo(i,:)=A1*exp(-1j*2*pi*(mu/2*(t+delt_t1(i)).^2)+-1j*2*pi*(fc)*(delt_t1(i)));
    echo(i,:)=awgn(echo(i,:),SNR);%%加噪声
    echo_fft(i,:)=(fft(echo(i,:)));%fftshift
    pc_result(i,:)=(ifft(echo_fft(i,:).*ht_fft));%ifftshift
    pc_result_fft(i,:)=(fft(pc_result(i,:)));%fftshift%快时间FFT
    pc_result(i,:)=(pc_result(i,:))/max(abs(pc_result(i,:)));%归一化
end
MTD = fft(pc_result,[],1);
[hang,lie]=find(max(max(abs(MTD)))==abs(MTD));
%%%%确定门限%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % MonteCarloPfa=1/Pfa*100;
% % vt = exp(-1j*2*pi*fd*(0:M-1).'*Tr);  %导向矢量
% % rhoR=0.95;
% % for i=1:M
% %     for j=1:M
% %         R(i,j)=rhoR^abs(i-j);
% %     end
% % end
% % iR=inv(R);
% % R_half=R^0.5;
% % for i=1:MonteCarloPfa
% %     i
% %     X=(randn(M,L)+1i*randn(M,L))/sqrt(2); % 产生方差为1的复高斯白噪声 % Rwhite1=1/snapshot1*X1*X1'; eig(Rwhite1); % round(mean(abs(eig(Rwhite1)))) == 1
% %     S=(R_half*X)*(R_half*X)'; % 有L个训练样本估计的杂波与噪声的协方差矩阵(Rhalf*X表示接收的L个训练数据)
% %     iS=inv(S);
% %     W=(randn(M,1)+1i*randn(M,1))/sqrt(2); % 1i == -i
% %     x=R_half*W;%+pp; % noise=(randn(N,1)+j*randn(N,1))/sqrt(2);  % 接收信号仅包括杂波和噪声
% %     Tamf(i)=abs(vt'*iS*x)^2/abs(vt'*iS*vt);     %%%%%% AMF
% %     MTD_t = fft(x);
% %     MTD_th(i) = abs(MTD_t(hang));
% % end
% % TAMF=sort(Tamf,'descend');
% % MTD_TH=sort(MTD_th,'descend');
% % Th_AMF=(TAMF(floor(MonteCarloPfa*Pfa-1))+TAMF(floor(MonteCarloPfa*Pfa)))/2;
% % Th_MTD=(MTD_TH(floor(MonteCarloPfa*Pfa-1))+MTD_TH(floor(MonteCarloPfa*Pfa)))/2;
%%MC开始%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR = -40:-10;
MC = 1000;
for i_SNR = 1:length(SNR)
    i_SNR
    count_MTD=0;
    count_AMF=0;
    for i_MC = 1:MC
        for i=1:pusle_num
            echo(i,:)=A1*exp(-1j*2*pi*(mu/2*(t+delt_t1(i)).^2)+-1j*2*pi*(fc)*(delt_t1(i)));
            echo(i,:)=awgn(echo(i,:),SNR(i_SNR));%%加噪声
            echo_fft(i,:)=(fft(echo(i,:)));%fftshift
            pc_result(i,:)=(ifft(echo_fft(i,:).*ht_fft));%ifftshift
            pc_result_fft(i,:)=(fft(pc_result(i,:)));%fftshift%快时间FFT
            pc_result(i,:)=(pc_result(i,:))/max(abs(pc_result(i,:)));%归一化
        end
        %%MTD检测
        MTD = fft(pc_result,[],1);
%         cankao_MTD = [MTD(:,1:lie-50),MTD(:,lie+50:end)];
% %         cankao_MTD = [MTD(hang-3,lie-3:lie+3),... 
% %                       MTD(hang-2,lie-3),MTD(hang-2,lie+3),...
% %                       MTD(hang-1,lie-3),MTD(hang-1,lie+3),...
% %                       MTD(hang,lie-3),MTD(hang,lie+3),...
% %                       MTD(hang+1,lie-3),MTD(hang+1,lie+3),...
% %                       MTD(hang+2,lie-3),MTD(hang+2,lie+3),...
% %                       MTD(hang+3,lie-3:lie+3)];
%         mean_MTD = mean(mean(abs(cankao_MTD)));
%         [num_hang,num_lie] = size(cankao_MTD);
%         alpha=(num_hang*num_lie)*(Pfa^(-1/(num_hang*num_lie))-1);%求alpha
%         Th_MTD=alpha*mean_MTD;%%求门限;
        if abs(MTD(hang,lie))>3.8871%Th_MTD
            count_MTD=count_MTD+1;
        end
        %%AMF检测
% %         %%AMFs
% %         cankao = [pc_result(:,1:L/2-5),pc_result(:,L/2+5:end)];
% %         R = (cankao*cankao')/length(cankao);
% %         inv_R = inv(R);
% %         fdt = 2*Vrt/lamda;
% %         AMF = zeros(M,L);
% %         for i_fd = 1:length(fdt)
% %             for j = 1:L
% %                 S = exp(-1j*2*pi*fdt(i_fd)*(0:M-1).'*Tr);
% %                 AMF(i_fd,j)=abs(abs(S'*inv_R*pc_result(:,j))^2/(S'*inv_R*S));
% %             end
% %         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cankao = [pc_result(:,1:L/2-5),pc_result(:,L/2+5:end)];
        R = (cankao*cankao')/length(cankao);
        inv_R = inv(R);
        AMF=abs(abs(S'*inv_R*pc_result(:,lie))^2/(S'*inv_R*S))/abs(pc_result(:,lie)'*inv_R*pc_result(:,lie));%ACE;%
        if AMF>fun_Th_ACE(245,64,Pfa)%0.0334%fun_Th_AMF(245,64,Pfa)%0.0328fun_Th_ACE
            count_AMF=count_AMF+1;
        end
    end
    Pd_MTD(i_SNR) = count_MTD/MC;
    Pd_AMF(i_SNR) = count_AMF/MC;
end
figure
plot(SNR,Pd_MTD)
hold on
plot(SNR,Pd_AMF,'r')


