clear all;  close all; clc
% 2010.5.31
% 刘维建
% 《空时自适应信号处理》，5.9节 
% 组合空时主通道的自适应处理（CMCAP）
N=16; % 空域通道数
K=18; % 脉冲数
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
Tphi0=kron(eye(K),at);   % (5.9.5)
for m=1:length(fd)
    b=exp(j*pi*(0:K-1)'*fd(m));
    Tf0=kron(b,eye(N));         % (5.9.4)
    T1=[Tphi0,Tf0(:,1:end-1)];  % (5.9.8)
    Ru=T1'*R*T1;
    vt=kron(b,at);
    gt=kron(tb.*b,ta.*at);
    gtBar=T1'*gt;
    vtBar=T1'*vt;
    W=inv(Ru)*gtBar;
    SINR_in=abs(vt'*vt)/(10^(CNR/10)+1);
	SINR_out=abs(W'*vtBar)^2/(W'*Ru*W);
    IF(m)=SINR_out/SINR_in;
    WOPT=invR*vt;
    SINR_OUT=abs(WOPT'*vt)^2/(WOPT'*R*WOPT);
    IFOPT(m)=SINR_OUT/SINR_in;
end
figure;plot(fd,10*log10(abs(IF)),'b-x',fd,10*log10(abs(IFOPT)),'r')
hold on
TheoreticalMax=10*log10((N*K)*(10^(CNR/10)+1));
plot([-1,1],[TheoreticalMax,TheoreticalMax],'m:')
xlabel('归一化doppler')
ylabel('IF（dB）')
title('改善因子')
h=legend('CMCAP','optimum','理论值')
% set(h,'Interpreter','none') % set(h,'Interpreter','none', 'Box', 'off')
% set(h, 'Box', 'off')

