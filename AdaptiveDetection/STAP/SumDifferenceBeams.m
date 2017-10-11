clear all;  close all; clc
% 2010.6.9
% 刘维建
% 《空时自适应信号处理》，5.9节 
% 和差波束法
N=16; % 空域通道数
K=18; % 脉冲数
q=7;
CNR=60;
R=fun_GenerateSimpleR(K,N,CNR);
invR=inv(R);
num1=101;
num2=51;
theta=linspace(-1,1,num1);
fd=linspace(-1,1,num2);
theta_s=0;
beta=1;
at=exp(j*pi*(0:N-1)'*theta_s);
ta=chebwin(N,30);
tb=chebwin(K,40);
Ts=[ones(N,1),(-1).^(0:N-1)'];

U=exp(-j*2*pi/K*(0:K-1)'*(0:K-1));
tf=chebwin(K,30);
F=diag(tf)*U;     %(227)
K1=K-q+1;
% for m=1:length(fd)-q+1    
%     b=exp(j*pi*(0:K-1)'*fd(m));
%     Tf=exp(j*pi*(0:K-1)'*fd(m+(0:q-1)));
for m=(q-1)/2+1:length(fd)-(q-1)/2    
    b=exp(j*pi*(0:K-1)'*fd(m));
    Tf=exp(j*pi*(0:K-1)'*fd(m+(-(q-1)/2:(q-1)/2)));      
%      Tf=toeplitz([F(1:K1,m);zeros(q-1,1)],[F(1,m),zeros(1,q-1)]);  F$A
    T=kron(Tf,Ts);
    Ru=T'*R*T;
    vt=kron(b,at);        %
    gt=kron(tb.*b,ta.*at); % gt=kron(tb.*bt,ta.*at);
    gtBar=T'*gt;
    W=inv(Ru)*gtBar;
    SINR_in=abs(vt'*vt)/(10^(CNR/10)+1);
	SINR_out=abs(W'*(T'*vt))^2/(W'*Ru*W);
    IF(m-(q-1)/2)=SINR_out/SINR_in;
    WOPT=invR*vt;
    SINR_OUT=abs(WOPT'*vt)^2/(WOPT'*R*WOPT);
    IFOPT(m-(q-1)/2)=SINR_OUT/SINR_in;
end
fdd=fd(1:length(fd)-q+1);
figure;plot(fdd,10*log10(abs(IF)),'b-x',fdd,10*log10(abs(IFOPT)),'r')
hold on
TheoreticalMax=10*log10((N*K)*(10^(CNR/10)+1));
plot([-1,1],[TheoreticalMax,TheoreticalMax],'m:')
xlabel('归一化doppler')
ylabel('IF（dB）')
title('改善因子')
h=legend('\Sigma-\Delta','optimum','理论值')