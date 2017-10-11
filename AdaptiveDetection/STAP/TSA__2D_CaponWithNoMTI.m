clear all;  close all; clc
% 2010.5.28
% 刘维建
% 《空时自适应信号处理》，5.4.1节 
% 未采用MTI的时空二维Capon法（TSA），即时空级联法中的先时域滤波后空域自适应波束形成
N=12; % 空域通道数
K=18; % 脉冲数
CNR=60
num=51;
theta=linspace(-1,1,num);
beta=1;
fd=beta*theta;
theta_s=0;
fd_s=0.6;
at=exp(j*pi*(0:N-1)'*theta_s);
bt=exp(j*pi*(0:K-1)'*fd_s);
ta=chebwin(N,30);
tb=chebwin(K,40);
td=chebwin(K-2,40);
vt=kron(bt,at); % gt=kron(tb.*bt,ta.*at);
Ac=(10^(CNR/10))^0.5; % 设噪声功率为1
Rc=zeros(K*N);
for i=1:length(theta)
    a=exp(j*pi*(0:N-1)'*theta(i));
    b=exp(j*pi*(0:K-1)'*beta*theta(i));
    v=kron(b,a);
    Rc=Rc+v*v';
end
Rc=Rc/sum(eig(Rc))*Ac^2; power=10*log10(abs(sum(eig(Rc))))
R=Rc+eye(K*N)/sum(eig(eye(K*N)));
a=[1,-2*exp(j*pi*beta*theta_s),exp(j*2*pi*beta*theta_s)];
% A=toeplitz([a(1);zeros(K-3,1)],[a,zeros(1,K-3)]); % A=kron(eye(K-2),a):这是错误的！
Wtk=zeros(K,K);
Wsk=zeros(N,K);
for k=1:K
    U=exp(j*2*pi/K*(0:K-1)'*k);
    Wtk(:,k)=tb.*U;
    T=kron(Wtk(:,k),eye(N));
    Rk=T'*R*T;
    Wsk(:,k)=inv(Rk)*(ta.*at);
    W=kron(Wtk(:,k),Wsk(:,k));
    SINRout(k)=abs(W'*vt)^2/(W'*R*W);
end
[MaxSINRout,index]=max(SINRout);
Wopt=kron(Wtk(:,index),Wsk(:,index));
% figure;plot(10*log10(abs(SINRout)),'b-x')
for i=1:length(theta)
    a=exp(j*pi*(0:N-1).'*theta(i));
    for jj=1:length(fd)        
        b=exp(j*pi*(0:K-1).'*fd(jj));
        steer=kron(b,a);
        res_opt(jj,i)=(Wopt'*steer).^2;   
    end
end
[Theta Fd]=meshgrid(theta,fd);
res_opt=res_opt/max(max(abs(res_opt)));
res_opt=10*log10(abs(res_opt));
res_opt(find(res_opt<-200))=-200;
figure;mesh(Theta,Fd,res_opt)   % 在方位向-0.9处形成深凹口！！！！！！！！
xlabel('方位');ylabel('Doppler');zlabel('功率 (dB)')
figure;contour(Theta,Fd,res_opt)
xlabel('方位');ylabel('Doppler');zlabel('功率 (dB)')
% -------------------------------------------------------
invR=inv(R);
for m=1:length(fd)
    bt=exp(j*pi*(0:K-1).'*fd(m)); 
    vt=kron(bt,at); 
    gt=kron(tb.*bt,ta.*at); 
    for k=1:K-2
        U=exp(j*2*pi/K*(0:K-1)'*k);
        Wtk(:,k)=tb.*U;
        T=kron(Wtk(:,k),eye(N));
        Rk=T'*R*T;
        Wsk(:,k)=inv(Rk)*(ta.*at);
        W=kron(Wtk(:,k),Wsk(:,k));
        SINRout(k)=abs(W'*vt)^2/(W'*R*W);
    end
    [MaxSINRout,index]=max(SINRout);
    Wopt=kron(Wtk(:,index),Wsk(:,index));
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
xlabel('归一化doppler')
ylabel('IF（dB）')
title('改善因子')
legend('TSA','optimum','理论值')


