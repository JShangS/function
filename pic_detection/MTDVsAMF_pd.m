clc
clear
close all
%%%MTD��AMF�ļ����ʱȽϣ�MTD���������CA-CFAR��⣬AMFֱ�Ӻ�����ֵ���бȽϵõ��������
fc=1000e6;%��Ƶ%2000
B=1e6;%����%4
Tao=127e-6;%����
Fs=2*B;%����Ƶ��
t=-Tao/2:1/Fs:Tao/2-1/Fs;%����ʱ��
L=length(t);
mu=B/Tao;%��Ƶ��
C=3e8;
R_max=C*Tao/2;
delt_R=C/(2*Fs);%%�������뵥Ԫ
R_start1=L/2*delt_R;%Ŀ��1��ʼλ��%170
lamda=C/fc;
PRF=10000;
Tr=1/PRF;
Vr_start1=100;%1��ʼ�ٶ�
fd = 2*(Vr_start1+3)/lamda;
pusle_num=64;%������
PV=PRF*lamda/4;
M=pusle_num;
dV = PV/M;
Vrt = -PV/2+dV:dV:PV/2;
echo=zeros(M,L);%�ز�
Pfa = 1e-6;    % �龯��
r0 = -log(Pfa); % AMF����<A new CFAR detection test for radar>������
S = exp(-1j*2*pi*fd*(0:M-1).'*Tr);  %����ʸ��
for i=1:pusle_num
    Vr1(i)=Vr_start1;
    fd1(i)=2*Vr1(i)/lamda;
    delt_t1(i)=2*(R_start1+Vr1(i)*Tr*(i-1))/C;%�ز��ӳ�
end
%%��ѹϵ��
    ht_t=exp(-1j*2*pi*(mu/2*(t).^2));
    ht=conj(fliplr(ht_t));
    ht_fft=fftshift(fft(ht));%fftshift
A1=1; 
SNR = 25
for i=1:pusle_num
    echo(i,:)=A1*exp(-1j*2*pi*(mu/2*(t+delt_t1(i)).^2)+-1j*2*pi*(fc)*(delt_t1(i)));
    echo(i,:)=awgn(echo(i,:),SNR);%%������
    echo_fft(i,:)=(fft(echo(i,:)));%fftshift
    pc_result(i,:)=(ifft(echo_fft(i,:).*ht_fft));%ifftshift
    pc_result_fft(i,:)=(fft(pc_result(i,:)));%fftshift%��ʱ��FFT
%     pc_result(i,:)=(pc_result(i,:))/max(abs(pc_result(i,:)));%��һ��
end
MTD = abs(fft(pc_result,[],1));
% MTD = awgn(MTD,SNR);
% MTD = MTD/max(max(MTD));
mesh(abs(MTD));
[hang,lie]=find(max(max(abs(MTD)))==abs(MTD));

%%MC��ʼ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR = -20:-10;
MC = 1000;
for i_SNR = 1:length(SNR)
    i_SNR
    count_MTD=0;
    count_AMF=0;
    for i_MC = 1:MC
        for i=1:pusle_num
            echo(i,:)=A1*exp(-1j*2*pi*(mu/2*(t+delt_t1(i)).^2)+-1j*2*pi*(fc)*(delt_t1(i)));
            echo(i,:)=awgn(echo(i,:),SNR(i_SNR));%%������
            echo_fft(i,:)=(fft(echo(i,:)));%fftshift
            pc_result(i,:)=(ifft(echo_fft(i,:).*ht_fft));%ifftshift
            pc_result_fft(i,:)=(fft(pc_result(i,:)));%fftshift%��ʱ��FFT
%             pc_result(i,:)=(pc_result(i,:))/max(abs(pc_result(i,:)));%��һ��
        end
        %%MTD���
        MTD = abs(fft(pc_result,[],1));
%         MTD = MTD/max(max(MTD));
%         MTD = awgn(MTD,SNR(i_SNR));
        cankao_MTD = [MTD(:,1:lie-50),MTD(:,lie+50:end)];
%         cankao_MTD = [MTD(hang-3,lie-3:lie+3),... 
%                       MTD(hang-2,lie-3),MTD(hang-2,lie+3),...
%                       MTD(hang-1,lie-3),MTD(hang-1,lie+3),...
%                       MTD(hang,lie-3),MTD(hang,lie+3),...
%                       MTD(hang+1,lie-3),MTD(hang+1,lie+3),...
%                       MTD(hang+2,lie-3),MTD(hang+2,lie+3),...
%                       MTD(hang+3,lie-3:lie+3)];
        mean_MTD = mean(mean(abs(cankao_MTD)));
        [num_hang,num_lie] = size(cankao_MTD);
        alpha=(num_hang*num_lie)*(Pfa^(-1/(num_hang*num_lie))-1);%��alpha
        Th_MTD=alpha*mean_MTD;%%������;
        if abs(MTD(hang,lie))>Th_MTD
            count_MTD=count_MTD+1;
        end
        %%AMF���
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
%         cankao = [pc_result(:,1:L/2-5),pc_result(:,L/2+5:end)];
%         K = length(cankao);
%         N = M;
%         R = (cankao*cankao')/length(cankao);
%         inv_R = inv(R);
%         AMF = abs(abs(S'*inv_R*pc_result(:,lie))^2/(S'*inv_R*S));
%         ACE = AMF/abs(pc_result(:,lie)'*inv_R*pc_result(:,lie));
%         if ACE> fun_Th_ACE(K,N,Pfa)   %AMF>r0%0.0328
%             count_AMF=count_AMF+1;
%         end
    end
    Pd_MTD(i_SNR) = count_MTD/MC;
%     Pd_AMF(i_SNR) = count_AMF/MC;
end
figure
plot(SNR,Pd_MTD,'b','LineWidth',2)
hold on
% plot(SNR,Pd_AMF,'r','LineWidth',2)







% PFA = 10e-3;
% MonteCarloPfa=1/PFA*100;
% vt = exp(-1j*2*pi*fd*(0:M-1).'*Tr);  %����ʸ��
% rhoR=0.95;
% for i=1:M
%     for j=1:M
%         R(i,j)=rhoR^abs(i-j);
%     end
% end
% iR=inv(R);
% R_half=R^0.5;
% for i=1:MonteCarloPfa
%     X=(randn(M,L)+1i*randn(M,L))/sqrt(2); % ��������Ϊ1�ĸ���˹������ % Rwhite1=1/snapshot1*X1*X1'; eig(Rwhite1); % round(mean(abs(eig(Rwhite1)))) == 1
%     S=(R_half*X)*(R_half*X)'; % ��L��ѵ���������Ƶ��Ӳ���������Э�������(Rhalf*X��ʾ���յ�L��ѵ������)
%     iS=inv(S);
%     W=(randn(M,1)+1i*randn(M,1))/sqrt(2); % 1i == -i
%     x=R_half*W;%+pp; % noise=(randn(N,1)+j*randn(N,1))/sqrt(2);  % �����źŽ������Ӳ�������
%     Tamf(i)=abs(vt'*iS*x)^2/abs(vt'*iS*vt);     %%%%%% AMF
% end
% TAMF=sort(Tamf,'descend');
% Th_AMF=(TAMF(floor(MonteCarloPfa*PFA-1))+TAMF(floor(MonteCarloPfa*PFA)))/2;