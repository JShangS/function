clear all;close all;clc
% 2010.3.17
% ��ά��
% ��Space-Time Adaptive Processing for Airborne Radar��,J.Ward, 5.3.1��
% Ԫ�ؿռ�������STAP

N=18;      %��Ԫ��
M=18;      %������
Num=500;
beta=1;
theta=(-1/2+1/Num:1/Num:1/2).';
theta0=linspace(-0.5,0.5,M);
fd=beta*theta;
CNR=40;
INR=40;
R =zeros(N*M); 
Rc=zeros(N*M); 
Ac=(10^(CNR/10))^0.5;
% Aj=(10^(INR/10))^0.5;
% theta_j=-0.3; %���ŵĹ�һ����λ
theta_s=theta0(10);  %Ŀ��Ĺ�һ����λ
fd_s=theta0(5);      %Ŀ��Ĺ�һ��������Ƶ��
at=exp(j*2*pi*(0:N-1).'*theta_s);
bt=exp(j*2*pi*(0:M-1).'*fd_s);
% ta=chebwin(N,30);
% tb=chebwin(M,30);
vt=kron(bt,at); %Ŀ���ʱ��������
% aj=exp(j*2*pi*(0:N-1).'*theta_j);
% Rj=kron(eye(M),Aj^2*aj*aj'); %��100��
for i=1:length(theta)
    a=exp(j*2*pi*(0:N-1).'*(theta(i)-theta_s));      %��31��
    b=exp(j*2*pi*(0:M-1).'*beta*(theta(i)-theta_s)); %��33��
    v=kron(b,a);
    Rc=Rc+v*v';
%     Rc=Rc+Ac^2*kron(b*b',a*a');
end
Rc0=Rc/max(max(abs(Rc)));
R=Ac^2*Rc0+eye(M*N);

% U=exp(-j*2*pi/M*(0:M-1)'*(0:M-1));
% U=exp(-j*2*pi*theta0'*2);
% % td=chebwin(M,30);
% F=conj(U);%*diag(td);     %(189)
% ta=chebwin(N,30);
% gt=ta.*at;          % P113�����ε�һ��
% vt=kron(bt,gt); %Ŀ���ʱ��������
wBar=zeros(N,M);
for p=1:M
    U=exp(-j*2*pi/M*p*(0:M-1)');
    T=kron(U,eye(N));
    Ru=T'*R*T;                 %(196)
    vts=T'*vt;
    wBar=inv(Ru)*vts;           %(195)
%     wm=kron(F(:,p),wBar(:,p));      %(202)
    SINR(p)=abs(wBar'*vts)^2/(wBar'*Ru*wBar);  %(185)  %�˴���vt��vt=kron(bt,at)������vt=kron(bt,gt)?
%     SINR_input=1/(1+Ac^2);
%     IF0(p)=SINR(p)/SINR_input;
end
SINR0=(abs(SINR));

% [Theta Fd]=meshgrid(theta,fd);
% res_opt=res_opt/max(max(abs(res_opt)));
% res_opt=10*log(abs(res_opt));
% IF=10*log(abs(IF));
% res_opt(find(res_opt<-250))=-250;
figure;plot(1:18,SINR0,'*-')
% xlabel('��λ');ylabel('Doppler');zlabel('���� (dB)')

% [SINRmax,index]=max(SINR);
% Wopt=kron(F(:,index),wBar(:,index));
% 
% for i=1:length(theta)
%     a=exp(j*2*pi*(0:N-1).'*theta(i));
%     for jj=1:length(fd)        
%         b=exp(j*2*pi*(0:M-1).'*fd(jj));
%         steer=kron(b,a);
% %         steer=steer/norm(steer);
%         res_opt(jj,i)=(steer'*Wopt).^2;        
% %         IF(i,jj)=steer'*inv_R*steer/(steer'*steer);
%         IF(jj,i)=(1+10^(CNR/10))*abs(Wopt'*steer)^2/(Wopt'*R*Wopt);        
%     end
% end
% 
% 
% for jj=1:length(fd)
%     b=exp(j*2*pi*(0:M-1).'*fd(jj))/sqrt(M);
%     steer=kron(b,at); 
% %     steer=steer/norm(steer);
%     IF0(jj)=(1+10^(CNR/10))*abs(Wopt'*steer)^2/(Wopt'*R*Wopt);        
% end
% IF0=10*log(abs(IF0));
% 
% [Theta Fd]=meshgrid(theta,fd);
% res_opt=res_opt/max(max(abs(res_opt)));
% res_opt=10*log(abs(res_opt));
% IF=10*log(abs(IF));
% res_opt(find(res_opt<-250))=-250;
% figure;mesh(Theta,Fd,res_opt)
% xlabel('��λ');ylabel('Doppler');zlabel('���� (dB)')
% figure;contour(Theta,Fd,res_opt)
% xlabel('��λ');ylabel('Doppler');zlabel('���� (dB)')
% figure;mesh(Theta,Fd,IF)
% xlabel('��λ');ylabel('Doppler');zlabel('���� (dB)')
% 
% figure;plot(fd,IF0)
% hold on
% plot([-0.5,0.5],[0,0],':')
% xlabel('��һ��doppler')
% ylabel('�������ӣ�dB��')
% title('\theta=\theta_s')
% 
% 
% 
