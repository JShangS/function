clear all;close all;clc
% 2010.5.22
% 刘维建
% 《空时自适应信号处理》，5.2.2节，先时域自适应滤波后空域自适应波束形成 

N=16; % 空域通道数
K=32; % 脉冲数
CNR=60;
num=101;
R=fun_GenerateSimpleR(K,N,num,CNR);
fd=linspace(-1,1,num);
theta_s=0.3;
at=exp(j*pi*(0:N-1)'*theta_s);
ta=chebwin(N,30);
tb=chebwin(K,30);
Beam=ta.*at;
Ac=(10^(CNR/10))^0.5; % 设噪声功率为1
vector_unit=[1;zeros(N-1,1)]; % vector_zero=zeros(N,1);
T2=kron(eye(K),vector_unit);
T1=T2.';
Ryt=T1*R*T2;
inv_Ryt=inv(Ryt);
invR=inv(R);
for m=1:length(fd)
    b=exp(j*pi*(0:K-1)'*fd(m));
    vt=kron(b,at); % gt=kron(tb.*bt,ta.*at);  
    F=tb.*b;
    wt=inv_Ryt*F;
    Rys=kron(wt,eye(N))'*R*kron(wt,eye(N));
    ws=inv(Rys)*Beam;
    Wopt=kron(wt,ws);
    SINR_in=abs(vt'*vt)/(10^(CNR/10)+1);
	SINR_out=abs(Wopt'*vt)^2/(Wopt'*R*Wopt);
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
h=legend('SA_TA','optimum','理论值')
set(h,'Interpreter','none') % set(h,'Interpreter','none', 'Box', 'off')
set(h, 'Box', 'off')
% ------------------------------------------------------------------------
        theta = linspace(-1,1,num);
        F=zeros(K);
        wt=zeros(K,K);
        ws=zeros(N,K);
        for k=1:K
            u=exp(-1j*2*pi/K*(0:K-1)'*k);
            F=tb.*u;
            wt(:,k)=inv_Ryt*F;
            Rys=kron(wt(:,k),eye(N))'*R*kron(wt(:,k),eye(N));
            ws(:,k)=inv(Rys)*Beam;
            w=kron(wt(:,k),ws(:,k));
            SINRout(k)=abs(w'*vt)^2/(w'*R*w);
        end
        [MaxSINRout,index]=max(SINRout);
        Wopt=kron(wt(:,index),ws(:,index));
        for i=1:length(theta)
            a=exp(j*pi*(0:N-1).'*theta(i));
            for jj=1:length(fd)        
                b=exp(j*pi*(0:K-1).'*fd(jj));
                steer=kron(b,a);
                res_opt(jj,i)=(steer'*Wopt).^2;   
            end
        end
        [Theta Fd]=meshgrid(theta,fd);
        res_opt=res_opt/max(max(abs(res_opt)));
        res_opt=10*log10(abs(res_opt)); % res_opt1(find(res_opt)<-200)=-200;
        figure;mesh(Theta,Fd,res_opt)
        xlabel('方位');ylabel('Doppler');zlabel('功率 (dB)')
        figure;contour(Theta,Fd,res_opt)
        xlabel('方位');ylabel('Doppler');zlabel('功率 (dB)')
