clear all;  close all; clc
% 2010.5.20
% 刘维建
% 《空时自适应信号处理》，5.2.1节 
% 常规空时级联（Conventional Space Time Cascade），即：先空域波束形成后时域滤波
N=16; % 空域通道数
K=32; % 脉冲数
Num=181;
CNR=60;
R=fun_GenerateSimpleR(K,N,Num,CNR);
num=101;
fd=linspace(-1,1,num);
theta_s=0;
at=exp(j*pi*(0:N-1)'*theta_s);
ta=chebwin(N,70); % ta=ones(N,1);
tb=chebwin(K,70); % tb=ones(K,1);
Beam=ta.*at;
invR=inv(R);
for m=1:length(fd)
    b=exp(j*pi*(0:K-1)'*fd(m));    
    vt=kron(b,at); %  gt=kron(tb.*bt,ta.*at);  
    NormVt=norm(vt);
    vtNorm=vt/norm(vt);
    F=tb.*b;
    Wopt=kron(F,Beam);
    SINR_in(m)=NormVt^2/(10^(CNR/10)+1);
	SINR_out(m)=abs(Wopt'*vt)^2/(Wopt'*R*Wopt);
%     SINR_out(m)=abs(gt'*invR*vt)^2/(gt'*invR*gt); %J. Ward报告中(113)式
    IF(m)=SINR_out(m)/SINR_in(m);
    WOPT=invR*vtNorm;
%     SINR_OUT(m)=NormVt^2*abs(WOPT'*vtNorm)^2/(WOPT'*R*WOPT);
    SINR_OUT(m)=NormVt^2*abs(vtNorm'*invR*vtNorm);
    IFOPT(m)=SINR_OUT(m)/SINR_in(m);
end
figure;plot(fd,10*log10(abs(IF)),'b-x',fd,10*log10(abs(IFOPT)),'r')
hold on
TheoreticalMax=10*log10((N*K)*(10^(CNR/10)+1));
plot([-1,1],[TheoreticalMax,TheoreticalMax],'m:')
xlabel('归一化doppler')
ylabel('IF（dB）')
title('改善因子')
% grid on
h=legend('CST','optimum','理论值')
% ------------------------------------------------------------------------
%     fd_s=0.3; % 不画频响图的话，fd_s用不着
%     bt=exp(j*pi*(0:K-1)'*fd_s);  % 不画频响图的话，bt用不着
%     vt=kron(bt,at); % gt=kron(tb.*bt,ta.*at);
%     F=zeros(K,K);
%     for k=1:K
%         u=exp(-j*2*pi/K*(0:K-1)'*k);
%         F(:,k)=tb.*u;
%         w=kron(F(:,k),Beam);
%         SINRout(k)=abs(w'*vt)^2/(w'*R*w);
%     end
%     [MaxSINRout,index]=max(SINRout);
%     Wopt=kron(F(:,index),Beam);
%     % figure;plot(1:K,abs(SINRout),'b-x');
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
%     res_opt=10*log10(abs(res_opt));
%     % res_opt1(find(res_opt)<-200)=-200;
%     figure;mesh(Theta,Fd,res_opt)   % 在方位向-0.9处形成深凹口！！！！！！！！
%     xlabel('方位');ylabel('Doppler');zlabel('功率 (dB)')
%     figure;contour(Theta,Fd,res_opt)
%     xlabel('方位');ylabel('Doppler');zlabel('功率 (dB)')


