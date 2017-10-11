clear all;  close all; clc
% 2010.5.20
% ��ά��
% ����ʱ����Ӧ�źŴ�����5.2.1��
% �ȿ��򳣹沨���γɺ�ʱ������Ӧ�˲���ST_ConventionalSpatialBeamforming_then_TemporalAdapt 
N=16; % ����ͨ����
K=32; % ������
CNR=60;
R=fun_GenerateSimpleR(K,N,CNR);

num=101;
fd=linspace(-1,1,num);
theta_s=0;
at=exp(j*pi*(0:N-1)'*theta_s);
ta=chebwin(N,60);
tb=ones(K,1); % tb=chebwin(K,60);
Beam=ta.*at;
Ac=(10^(CNR/10))^0.5; % ����������Ϊ1
Ry=kron(eye(K),Beam)'*R*kron(eye(K),Beam);
inv_Ry=inv(Ry);
invR=inv(R);
for m=1:length(fd)
    b=exp(j*pi*(0:K-1)'*fd(m));
    vt=kron(b,at);  % gt=kron(tb.*bt,ta.*at);
    F=tb.*b; % �������Ȩ��ʱ�򲻼�Ȩ��Ҫ�ȿ���ʱ�����Ȩ(��:tb=ones(K,1);)Ч��Ҫ�ã��� 
    wt=inv_Ry*F;
    Wopt=kron(wt,Beam);
    SINR_in(m)=abs(vt'*vt)/(10^(CNR/10)+1);
	SINR_out(m)=abs(Wopt'*vt)^2/(Wopt'*R*Wopt);
    IF(m)=SINR_out(m)/SINR_in(m);
    WOPT=invR*vt;
    SINR_OUT(m)=abs(WOPT'*vt)^2/(WOPT'*R*WOPT);
    IFOPT(m)=SINR_OUT(m)/SINR_in(m);
end
figure;plot(fd,10*log10(abs(IF)),'b-x',fd,10*log10(abs(IFOPT)),'r')
hold on
TheoreticalMax=10*log10((N*K)*(10^(CNR/10)+1));
plot([-1,1],[TheoreticalMax,TheoreticalMax],'m:')
xlabel('��һ��doppler')
ylabel('IF��dB��')
title('��������')
h=legend('SB_TA','optimum','����ֵ')
set(h,'Interpreter','none') % set(h,'Interpreter','none', 'Box', 'off')
set(h, 'Box', 'off')
% ------------------------------------------------------------------------
%     fd_s=0.3; % ����Ƶ��ͼ�Ļ���fd_s�ò���
%     bt=exp(j*pi*(0:K-1)'*fd_s);  % ����Ƶ��ͼ�Ļ���bt�ò���
%     vt=kron(bt,at); % gt=kron(tb.*bt,ta.*at);
%     wt=zeros(K);
%     for k=1:K
%         u=exp(-j*2*pi/K*(0:K-1)'*k);
%         F=tb.*u;
%         wt(:,k)=inv_Ry*F;
%         w=kron(wt(:,k),Beam);
%         SINRout(k)=abs(w'*vt)^2/(w'*R*w);
%     end
%     [MaxSINRout,index]=max(SINRout);
%     Wopt=kron(wt(:,index),Beam);
%     for i=1:length(theta)
%         a=exp(j*pi*(0:N-1).'*theta(i));
%         for jj=1:length(fd)        
%             b=exp(j*pi*(0:K-1).'*fd(jj));
%             steer=kron(b,a);
%             res_opt(jj,i)=(steer'*Wopt).^2;   
%         end
%     end
%     [Theta Fd]=meshgrid(theta,fd);
%     res_opt=res_opt/max(max(abs(res_opt)));
%     res_opt=10*log10(abs(res_opt)); % res_opt1(find(res_opt)<-200)=-200;
%     figure;mesh(Theta,Fd,res_opt) % �ڷ�λ��-0.9���γ���ڣ���������������
%     xlabel('��λ');ylabel('Doppler');zlabel('���� (dB)')
%     figure;contour(Theta,Fd,res_opt)
%     xlabel('��λ');ylabel('Doppler');zlabel('���� (dB)')
