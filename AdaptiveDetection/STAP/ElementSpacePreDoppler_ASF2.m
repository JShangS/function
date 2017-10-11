clear all; close all;clc
% 2010.3.16
% ��ά��
% ��Space-Time Adaptive Processing for Airborne Radar��,J.Ward, 5.2��
% Ԫ�ؿռ�Ԥ������STAP��A$F��
% δ�ɹ�

N=18;      %��Ԫ��
M=18;      %������
K=3;       %ÿ����CPI����K������
M0=M-K+1;  %��CPI�ĸ���
Num=100;
beta=1;
theta=(-1/2+1/Num:1/Num:1/2).';
fd=beta*theta;
CNR=40;
INR=40;
R =zeros(N*M); 
Rc=zeros(N*M); 
Ac=(10^(CNR/10))^0.5;
Aj=(10^(INR/10))^0.5;
theta_j=-0.3; %���ŵĹ�һ����λ
theta_s=0.1;  %Ŀ��Ĺ�һ����λ
fd_s=0.25;    %Ŀ��Ĺ�һ��������Ƶ��
at=exp(j*2*pi*(0:N-1).'*theta_s);
bt=exp(j*2*pi*(0:M-1).'*fd_s);
bp=exp(j*2*pi*(0:K-1).'*fd_s);
ta=chebwin(N,30);
tb=chebwin(M,30);
vt=kron(bt,at); %Ŀ���ʱ��������
% vt_weight=kron(tb,ta).*vt;  %(108)
tgb=[1;2;1];
ta=chebwin(N,30);
% gt=kron(tgb.*bp,ta.*at);       %(172)
gt=kron(bt,ta.*at);       %(172)
aj=exp(j*2*pi*(0:N-1).'*theta_j);
Rj=kron(eye(M),Aj^2*aj*aj'); %��100��
for i=1:length(theta)
    a=exp(j*2*pi*(0:N-1).'*theta(i));      %��31��
    b=exp(j*2*pi*(0:M-1).'*beta*theta(i)); %��33��
    v=kron(b,a);
    v=v/norm(v);
    Rc=Rc+Ac^2*v*v';
%     Rc=Rc+Ac^2*kron(b*b',a*a');
end
R=Rc+Rj+eye(M*N); % R=Rc+eye(M*N);
W=zeros(M*N,M0);
Jp=zeros(M,K);
W=zeros(K*N,M0);
for p=1:M0
%     Jp=[zeros(p-1,K);eye(K);zeros(M-K-p+1,K)];     %��175��
    Jp=zeros(M,K);
    Jp((1:K)+(p-1),:)=eye(K);
    T=kron(Jp,eye(N));
    Ru=T'*R*T;
    Jptaper=zeros(M,K);
    Jptaper((1:K)+(p-1),:)=diag([1,2,1]).*eye(K); % ����ʽ��Ȩ
    Ttaper=kron(Jptaper,eye(N));    
    gtBar=Ttaper'*gt;     
%     gt=kron(bt,ta.*at); gtBar=T'*gt;
    W(:,p)=inv(Ru)*gtBar;           %(168)
    SINR(p)=abs( W(:,p)'*gtBar)^2/( W(:,p)'*Ru* W(:,p));  %(185)
end
[SINRmax,index]=max(SINR);
Wopt=W(:,index);
Jp=zeros(M,K);
Jp((1:K)+(index-1),:)=eye(K);
T=kron(Jp,eye(N));
Wopt=T*Wopt;
for i=1:length(theta)
    a=exp(j*2*pi*(0:N-1).'*theta(i));
    for jj=1:length(fd)        
        b=exp(j*2*pi*(0:M-1).'*fd(jj));
        steer=kron(b,a);
        res_opt(jj,i)=(steer'*Wopt).^2;        
    end
end
[Theta Fd]=meshgrid(theta,fd);
res_opt=res_opt/max(max(abs(res_opt)));
res_opt=10*log(abs(res_opt));
res_opt(find(res_opt<-250))=-250;
figure;mesh(Theta,Fd,res_opt)
xlabel('��λ');ylabel('Doppler');zlabel('���� (dB)')
figure;contour(Theta,Fd,res_opt)
xlabel('��λ');ylabel('Doppler');zlabel('���� (dB)')
% -------------------------------------------------------
for n=1:length(fd)
    bt=exp(j*2*pi*(0:M-1).'*fd(n)); 
    vt=kron(bt,at);
%     vt=kron(tb.*bt,ta.*at);
    bp=exp(j*2*pi*(0:K-1).'*fd(n));
%     gt=kron(bp,at);       %(172)
%     gt=kron(tgb.*bp,at);     %(172)
    gt=kron(tgb.*bp,ta.*at);       %(172)
    for p=1:M0
        Jp=[zeros(p-1,K);eye(K);zeros(M-K-p+1,K)];     %��175��
        T=kron(Jp,eye(N));
        Ru=T'*R*T;
        w=inv(Ru)*gt;           %(168)
        W(1+(p-1)*N:K*N+(p-1)*N,p)=w;   %(180)
    end
    for p=1:M0
        wm=W*F(:,p);    %(184)
        SINR(p)=abs(wm'*vt)^2/(wm'*R*wm);  %(185), SINR�ǳ�С������������
    end
    [SINRmax,index]=max(SINR);
    Wopt=W*F(:,index);
	SINR_in(n)=abs(vt'*vt)/(10^(CNR/10)+1);
	SINR_out(n)=abs(Wopt'*vt)^2/(Wopt'*R*Wopt);
    IF(n)=SINR_out(n)/SINR_in(n);    
end
figure;plot(fd,10*log(abs(IF)))
xlabel('��һ��doppler')
ylabel('IF��dB��')
title('��������')
% figure;plot(fd,SINR_in)
% figure;plot(fd,SINR_out)