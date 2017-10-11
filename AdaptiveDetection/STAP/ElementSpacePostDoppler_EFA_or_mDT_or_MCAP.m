clear all;close all;clc
% 2010.3.17
% ��ά��
% ��Space-Time Adaptive Processing for Airborne Radar��,J.Ward, 5.3.2��
% Ԫ�ؿռ�������STAP��EFA����mDT

N=18;      %��Ԫ��
M=18;      %������
K=5;       %ÿ����CPI����K������
P=(K-1)/2; %(231)�Ϸ�����     
M0=M-K+1;  %��CPI�ĸ���
Num=100;
beta=1;
theta=(-1/2+1/Num:1/Num:1/2).';
fd=beta*theta;
CNR=40;
JNR=40;
R =zeros(N*M); 
Rc=zeros(N*M); 
Ac=(10^(CNR/10))^0.5;
Aj=(10^(JNR/10))^0.5;
theta_j=-0.3; %���ŵĹ�һ����λ
theta_s=0.1;  %Ŀ��Ĺ�һ����λ
fd_s=0.25;    %Ŀ��Ĺ�һ��������Ƶ��
at=exp(j*2*pi*(0:N-1).'*theta_s);
bt=exp(j*2*pi*(0:M-1).'*fd_s);
bp=exp(j*2*pi*(0:K-1).'*fd_s);
% ta=chebwin(N,30);
% tb=chebwin(M,30);
vt=kron(bt,at); %Ŀ���ʱ��������
aj=exp(j*2*pi*(0:N-1).'*theta_j);
aj=aj/norm(aj);
Rj=kron(eye(M),Aj^2*aj*aj'); %��100��
for i=1:length(theta)
    a=exp(j*2*pi*(0:N-1).'*theta(i));      %��31��
    b=exp(j*2*pi*(0:M-1).'*beta*theta(i)); %��33��
    v=kron(b,a);
    v=v/norm(v);
    Rc=Rc+Ac^2*v*v';
%     Rc=Rc+Ac^2*kron(b*b',a*a');
end
Rc=Rc/length(theta); power1=10*log10(abs(sum(eig(Rc))))
R=Rc+Rj+eye(M*N)/(M*N); % R=Rc+eye(K*N)/sum(eig(eye(K*N)));

U=exp(-j*2*pi/M*(0:M-1)'*(0:M-1));
td=chebwin(M,30);
F=diag(td)*conj(U);     %(230)

Ww=zeros(M*N,M0);
for m=1+P:M-P
    Fm=F(:,m-P:m+P);    %��231��
    T=kron(Fm,eye(N));
    Ru=T'*R*T;
    gt=T'*vt;                       %(207)
    w=inv(Ru)*gt;                   %(205)
    Ww(:,m)=T*w;                    %(209)
    SINR(m)=abs(Ww(:,m)'*vt)^2/(Ww(:,m)'*R*Ww(:,m));  %(185)
end
[SINRmax,index]=max(SINR);
Wopt=Ww(:,index);

for i=1:length(theta)
    a=exp(j*2*pi*(0:N-1).'*theta(i));
    for jj=1:length(fd)        
        b=exp(j*2*pi*(0:M-1).'*fd(jj));
        steer=kron(b,a);
        steer=steer/norm(steer);
        res_opt(jj,i)=(steer'*Wopt).^2;        
    end
end
[Theta Fd]=meshgrid(theta,fd);
res_opt=res_opt/max(max(abs(res_opt)));
res_opt=10*log10(abs(res_opt));
% res_opt(find(res_opt<-250))=-250;
figure;mesh(Theta,Fd,res_opt)
xlabel('��λ');ylabel('Doppler');zlabel('���� (dB)')
figure;contour(Theta,Fd,res_opt)
xlabel('��λ');ylabel('Doppler');zlabel('���� (dB)')

% -------------------------------------------------------
Ww=zeros(M*N,M0);
invR=inv(R);
SINR_in=abs(vt'*vt)/(10^(CNR/10)+10^(JNR/10)+1);
for n=1:length(fd)
    bt=exp(j*2*pi*(0:M-1).'*fd(n)); 
    vt=kron(bt,at);
    for m=1+P:M-P
        Fm=F(:,m-P:m+P);    %��231��
        T=kron(Fm,eye(N));
        Ru=T'*R*T;
        gt=T'*vt;                       %(207)
        w=inv(Ru)*gt;                   %(205)
        Ww(:,m)=T*w;                    %(209)
        SINR(m)=abs(Ww(:,m)'*vt)^2/(Ww(:,m)'*R*Ww(:,m));  %(185)
    end
    [SINRmax,index]=max(SINR);
    Wopt=Ww(:,index);
	SINR_out(n)=abs(Wopt'*vt)^2/(Wopt'*R*Wopt);
    IF(n)=SINR_out(n)/SINR_in;    
    WOPT=invR*vt;
    SINR_OUT(n)=abs(WOPT'*vt)^2/(WOPT'*R*WOPT);
    IFOPT(n)=SINR_OUT(n)/SINR_in;
end
figure;plot(fd,10*log10(abs(IF)),'b-x',fd,10*log10(abs(IFOPT)),'r-*')
hold on
TheoreticalMax=10*log10((M*N)*(10^(CNR/10)+10^(JNR/10)+1));
plot([-0.5,0.5],[TheoreticalMax,TheoreticalMax],'k:')
xlabel('��һ��doppler')
ylabel('IF��dB��')
title('��������')
h = legend('EFA','optimum','����ֵ',3)
set(h,'Interpreter','none') % set(h,'Interpreter','none', 'Box', 'off')
set(h, 'Box', 'off')
