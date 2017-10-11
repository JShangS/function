clear all;close all;clc
% 2010.3.17
% 刘维建
% 《Space-Time Adaptive Processing for Airborne Radar》,J.Ward, 5.3.1节
% 元素空间后多普勒STAP

N=18;      %阵元数
M=18;      %脉冲数
K=3;       %每个子CPI包括K个脉冲
% M0=M-K+1;  %子CPI的个数
Num=128;
beta=1;
theta=(-1/2+1/Num:1/Num:1/2).';
fd=beta*theta;
CNR=40;
JNR=40;
R =zeros(N*M); 
Rc=zeros(N*M); 
Ac=(10^(CNR/10))^0.5;
Aj=(10^(JNR/10))^0.5;
theta_j=-0.3; %干扰的归一化方位
theta_s=0.1;  %目标的归一化方位
fd_s=0.25;    %目标的归一化多普勒频率
at=exp(j*2*pi*(0:N-1).'*theta_s);
bt=exp(j*2*pi*(0:M-1).'*fd_s);
bp=exp(j*2*pi*(0:K-1).'*fd_s);
% ta=chebwin(N,30);
% tb=chebwin(M,30);
vt=kron(bt,at); %目标空时导向向量
aj=exp(j*2*pi*(0:N-1).'*theta_j);
aj=aj/norm(aj);
Rj=kron(eye(M),Aj^2*aj*aj'); %（100）
abs(10*log10(sum(eig(Aj^2*aj*aj'))))==JNR
for i=1:length(theta)
    a=exp(j*2*pi*(0:N-1).'*theta(i));      %（31）
    b=exp(j*2*pi*(0:M-1).'*beta*theta(i)); %（33）
    v=kron(b,a);
    v=v/norm(v);
    Rc=Rc+Ac^2*v*v';  % Rc=Rc+Ac^2*kron(b*b',a*a');
end
Rc=Rc/length(theta); power1=10*log10(abs(sum(eig(Rc))))
R=Rc+Rj+eye(M*N)/(M*N); % R=Rc+eye(K*N)/sum(eig(eye(K*N)));
U=exp(-j*2*pi/N*(0:M-1)'*(0:M-1));
td=chebwin(M,30);
F=diag(td)*conj(U);     %(189)
ta=chebwin(N,30);
% gt=ta.*at;
wBar=zeros(N,M);
for p=1:M
    Jp=[zeros(p-1,K);eye(K);zeros(M-K-p+1,K)];     %（175）
    T=kron(F(:,p),eye(N));
    Ru=T'*R*T;              %(196)
    vtBar=T'*vt;
 %     vtBar=F(:,p)'*bt*at;       %(192)~(194)    
    wBar(:,p)=inv(Ru)*vtBar;           %(195)
%     wm=kron(F(:,p),w);      %(202)
%     SINR(p)=abs(wm'*vt)^2/(wm'*R*wm);  %(185)
    SINR(p)=abs(wBar(:,p)'*vtBar)^2/(wBar(:,p)'*Ru*wBar(:,p));        
%     U=exp(-j*2*pi/M*p*(0:M-1)');
%     T=kron(U,eye(N));
%     Ru=T'*R*T;                 %(196)
%     vts=T'*vt;
%     wBar=inv(Ru)*vts;           %(195)
% %     wm=kron(F(:,p),wBar(:,p));      %(202)
%     SINR(p)=abs(wBar'*vts)^2/(wBar'*Ru*wBar);  %(185)  %此处的vt是vt=kron(bt,at)，还是vt=kron(bt,gt)?    
end
figure;plot(abs(SINR),'x-')
[SINRmax,index]=max(SINR);
Wopt=kron(F(:,index),wBar(:,index));

for i=1:length(theta)
    a=exp(j*2*pi*(0:N-1).'*theta(i));
    for jj=1:length(fd)        
        b=exp(j*2*pi*(0:M-1).'*fd(jj));
        steer=kron(b,a);    %   steer=steer/norm(steer);
        res_opt(jj,i)=(steer'*Wopt).^2;        
    end
end
[Theta Fd]=meshgrid(theta,fd);
res_opt=res_opt/max(max(abs(res_opt)));
res_opt=10*log10(abs(res_opt));
% res_opt(find(res_opt<-250))=-250;
figure;mesh(Theta,Fd,res_opt)
xlabel('方位');ylabel('Doppler');zlabel('功率 (dB)')
figure;contour(Theta,Fd,res_opt)
xlabel('方位');ylabel('Doppler');zlabel('功率 (dB)')
% -------------------------------------------------------
invR=inv(R);
for n=1:length(fd)
    b=exp(j*2*pi*(0:M-1).'*fd(n));
    vt=kron(b,at);     
    for p=1:M
        Jp=[zeros(p-1,K);eye(K);zeros(M-K-p+1,K)];     %（175）
        T=kron(F(:,p),eye(N));
        Ru=T'*R*T;              %(196)
        gt=T'*vt;
        wBar(:,p)=inv(Ru)*gt;           %(195)
        SINR(p)=abs(wBar(:,p)'*gt)^2/(wBar(:,p)'*Ru*wBar(:,p));        
    end
    [SINRmax,index]=max(SINR);
    Wopt=kron(F(:,index),wBar(:,index)); 
	SINR_in(n)=abs(vt'*vt)/(10^(CNR/10)+10^(JNR/10)+1);
	SINR_out(n)=abs(Wopt'*vt)^2/(Wopt'*R*Wopt);
    IF(n)=SINR_out(n)/SINR_in(n);    
    WOPT=invR*vt;
    SINR_OUT(n)=abs(WOPT'*vt)^2/(WOPT'*R*WOPT);
    IFOPT(n)=SINR_OUT(n)/SINR_in(n);
end
figure;plot(fd,10*log10(abs(IF)),'b-x',fd,10*log10(abs(IFOPT)),'r')
hold on
TheoreticalMax=10*log10((M*N)*(10^(CNR/10)+10^(JNR/10)+1));
plot([-0.5,0.5],[TheoreticalMax,TheoreticalMax],'m:')
xlabel('归一化doppler')
ylabel('IF（dB）')
title('改善因子')
h = legend('1DT','optimum','理论值',3)
set(h,'Interpreter','none') % set(h,'Interpreter','none', 'Box', 'off')
set(h, 'Box', 'off')

