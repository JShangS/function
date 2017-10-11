clear all;  close all; clc
% 2010.5.31
% 刘维建
% 《空时自适应信号处理》，5.10节 
% 广义相邻多波束（GMB）
N=16; % 空域通道数
K=18; % 脉冲数
CNR=60;
p=5; % 空域波束数     当p>7性能波动。
q=6; % 多普勒通道数   当p>7性能波动。
num1=101;
num2=51;
theta=linspace(-1,1,num1);
fd=linspace(-1,1,num2);
theta_s=0;
beta=1;
at=exp(j*pi*(0:N-1)'*theta_s);
ta=chebwin(N,30);
tb=chebwin(K,40);
Ac=(10^(CNR/10))^0.5; % 设噪声功率为1
Rc=zeros(K*N);
for i=1:length(theta)
    ac=exp(j*pi*(0:N-1)'*theta(i));
    bc=exp(j*pi*(0:K-1)'*beta*theta(i));
    v=kron(bc,ac);
    Rc=Rc+v*v';
end
Rc=Rc/sum(eig(Rc))*Ac^2; power=10*log10(abs(sum(eig(Rc))));
R=Rc+eye(K*N)/sum(eig(eye(K*N)));
invR=inv(R);
Tphi0=kron(eye(K),at);   % (5.9.4)
delta_theta=theta(2)-theta(1);
if mod(p,2)==0
    INDEX=-p/2+1:p/2;
else 
    INDEX=-(p-1)/2:(p-1)/2;
end
[Min,IndexOfTheta_s]=min(abs(theta_s-theta));
% Ts=exp(j*pi*(0:N-1)'*(theta_s+(-floor(p/2):floor(p/2))*delta_theta));
Ts=exp(j*pi*(0:N-1)'*theta(IndexOfTheta_s+INDEX));
for m=1:length(fd)-q+1    
    b=exp(j*pi*(0:K-1)'*fd(m));
    Tf=exp(j*pi*(0:K-1)'*fd(m+(0:q-1)));
    T=kron(Tf,Ts);
    Ru=T'*R*T;
    vt=kron(b,at);        %
    gt=kron(tb.*b,ta.*at); % gt=kron(tb.*bt,ta.*at);
    gtBar=T'*gt;
    W=inv(Ru)*gtBar;
    SINR_in=abs(vt'*vt)/(10^(CNR/10)+1);
	SINR_out=abs(W'*(T'*vt))^2/(W'*Ru*W);
    IF(m)=SINR_out/SINR_in;
    WOPT=invR*vt;
    SINR_OUT=abs(WOPT'*vt)^2/(WOPT'*R*WOPT);
    IFOPT(m)=SINR_OUT/SINR_in;
end
fdd=fd(1:length(fd)-q+1);
figure;plot(fdd,10*log10(abs(IF)),'b-x',fdd,10*log10(abs(IFOPT)),'r')
hold on
TheoreticalMax=10*log10((N*K)*(10^(CNR/10)+1));
plot([-1,1],[TheoreticalMax,TheoreticalMax],'m:')
xlabel('归一化doppler')
ylabel('IF（dB）')
title('改善因子')
h=legend('CMCAP','optimum','理论值')
% set(h,'Interpreter','none') % set(h,'Interpreter','none', 'Box', 'off')
% set(h, 'Box', 'off')

