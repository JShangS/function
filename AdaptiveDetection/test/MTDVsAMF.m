clc
clear
close all

%%%MTD��AMF�Ļ��۽����MTD���������CA-CFAR��⣬AMFֱ�Ӻ�����ֵ���бȽϵõ��������
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
Vr_start1=290;%1��ʼ�ٶ�
pusle_num=64;%������
PV=PRF*lamda/4;
M=pusle_num;
Pfa = 1e-6;     %�龯��
r0 = -log(Pfa); %����
echo=zeros(M,L);%�ز�
for i=1:pusle_num
    Vr1(i)=Vr_start1;
    fd1(i)=2*Vr1(i)/lamda;
    delt_t1(i)=2*(R_start1+Vr1(i)*Tr*(i-1))/C;%�ز��ӳ�
end
%%��ѹϵ��
    ht_t=exp(-1j*2*pi*(mu/2*(t).^2));
    ht=conj(fliplr(ht_t));
    ht_fft=fftshift(fft(ht));%fftshift
    SNR=-15;
A1=1; %%1.8%1.45
for i=1:pusle_num
    echo(i,:)=A1*exp(-1j*2*pi*(mu/2*(t+delt_t1(i)).^2)+-1j*2*pi*(fc)*(delt_t1(i)));
    echo(i,:)=awgn(echo(i,:),SNR);%%������
    echo_fft(i,:)=(fft(echo(i,:)));%fftshift
    pc_result(i,:)=(ifft(echo_fft(i,:).*ht_fft));%ifftshift
    pc_result_fft(i,:)=(fft(pc_result(i,:)));%fftshift%��ʱ��FFT
end
% % ��ѹ��ԭʼ�ز�
figure(1), 
mesh(abs(pc_result))
% imagesc (abs(pc_result))
xlabel('���뵥Ԫ')
ylabel('������')
%%MTD
dV = PV/M;
Vrt = -PV/2+dV:dV:PV/2;
tic
MTD = fft(pc_result,[],1);
toc
[X,Y] = meshgrid(1:L,Vrt);
figure(2)
mesh(X,Y,abs(MTD))
xlabel('���뵥Ԫ')
ylabel('�ٶ�m/s')
zlabel('����')
cankao_MTD = [MTD(:,1:L/2-50),MTD(:,L/2+50:end)];
mean_MTD = mean(mean(abs(cankao_MTD)));
[num_hang,num_lie] = size(cankao_MTD);
alpha=(num_hang*num_lie)*(Pfa^(-1/(num_hang*num_lie))-1);%��alpha
Th_MTD=alpha*mean_MTD;%%������;
bool_MTD = abs(MTD)>Th_MTD;
MTD_th = abs(MTD).*bool_MTD;
figure(3)
mesh(X,Y,MTD_th)
% imagesc(1:L,Vrt,abs(MTD))
%%AMFs
cankao = [pc_result(:,1:L/2-5),pc_result(:,L/2+5:end)];
R = (cankao*cankao')/length(cankao);
inv_R = inv(R);
fdt = 2*Vrt/lamda;
AMF = zeros(M,L);
tic
for i_fd = 1:length(fdt)
    for j = 1:L
        S = exp(-1j*2*pi*fdt(i_fd)*(0:M-1).'*Tr);
        AMF(i_fd,j)=abs(abs(S'*inv_R*pc_result(:,j))^2/(S'*inv_R*S));
        ACE(i_fd,j)=AMF(i_fd,j)/abs(pc_result(:,L/2)'*inv_R*pc_result(:,L/2));
    end
end
toc
figure(4)
mesh(X,Y,AMF)
xlabel('���뵥Ԫ')
ylabel('�ٶ�m/s')
zlabel('����')
N = M;
K = length(cankao);
bool_ACE = ACE>fun_Th_ACE(K,N,Pfa);
ACE_th = ACE.*bool_ACE;
figure(5)
mesh(X,Y,ACE_th)
% imagesc(1:L,Vrt,AMF)
% figure()
% plot(abs(MTD(:,L/2)));
% figure()
% plot(abs(AMF(:,L/2)));

