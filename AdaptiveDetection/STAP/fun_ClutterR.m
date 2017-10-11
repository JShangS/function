function R=fun_ClutterR(V,H,N,M,K,lambda,d,fr,phi0,theta0,Inn,Imm,Ikk,CNR,senser_error,channel_error,ea,ep,crab_theta,Br)
% fun_ClutterR is to generate covariance of the clutter.

V=140;                                            % 载机速度
H=8000;                                           % 载机高度
N=16;                                             % 每行有N个阵元
M=8;                                              % 每列有M个阵元
K=16;                                             % 时域采样数(时域脉冲数)
lambda=0.23;                                      % 波长
d=1/2*lambda;                                     % 阵元间距与波长之比
c=3e+8;                                           % 光速
fr=2434.8;                                        % PRF
Rmax=400000;                                      % 雷达最大作用距离
phi0=0;                                           % 波束中心高低角
theta0=pi/2;                                      % 波束中心方位角
Inn=30;                                           % 行向加权系数中衰减dB数
Imm=30;                                           % 列向加权系数中衰减dB数
Ikk=30;                                           % 时域加权系数中衰减dB数
CNR=50;                                           % 杂噪比

senser_error=0.0025;                              % 阵元误差百分比
channel_error=0.003;                              % 通道误差百分比
ea=0.0025*randn(M,N);                             % 阵元幅度误差
ep=0.0005*randn(M,N);                             % 阵元相位误差
crab_theta=0.1*pi/180;                            % 飞机偏航角
Br=0.005;                                         % 杂波起伏度
psi=0:pi/180:pi;                                  % 锥角取值范围
phi=0:pi/180:pi/2;                                % 高低角取值范围
theta=0:pi/180:pi;                                % 方位角取值范围
theta_step=theta(2)-theta(1);
caa=(channel_error*randn(N,1));ca=caa*caa.';                 %通道间幅度误差(列间幅度误差)
ca=ca-diag(diag(ca));
cpp=(channel_error*pi*randn(N,1));cp=cpp*cpp.';              %通道间相位误差(列间相位误差)
cp=cp-diag(diag(cp));
C=(1+ca).*exp(j*cp);    % (3.3.21)

Ru=c/(2*fr);
if  Ru >= H
     L=floor(Rmax/Ru)+1;
else L=floor(Rmax/Ru);
end

In=chebwin(N,Inn);
Im=chebwin(M,Imm);
Ik=chebwin(K,Ikk);
Rc=zeros(N*K,N*K);
Rc_el=zeros(N*K,N*K);
for el=0:L-1
    phi_el=asin(H/(H+el*Ru));
    Rel=1.5*H+el*Ru;   % 距离环的起始位置的取值范围为：[H,H+Ru]，此处设为：1.5H，
    for nn=1:length(theta)
        F(nn)=sum(In.*exp(j*2*pi*d/lambda*(0:N-1)'*(cos(theta(nn))*cos(phi_el)-cos(theta0)*cos(phi0))))*...
          sum(Im.*exp(j*2*pi*d/lambda*(0:M-1)'*(sin(phi_el)-sin(phi0))));     %(3.2.3)
    end
    omega_s=2*pi*d/lambda*cos(theta+crab_theta)*cos(phi_el);                         %(3.3.7a)
    omega_d=4*pi*V/(lambda*fr)*cos(theta+crab_theta)*cos(phi_el);                    %(3.3.24)
    for n1=1:N
        G(n1)=sum(Im.*(1+ea(:,n1)).*exp(j*ep(:,n1)).*exp(j*2*pi*d/lambda*((0:M-1)')*(sin(phi_el)-sin(phi0)))); %(3.2.20) 
        for n2=1:N
            G(n2)=sum(Im.*(1+ea(:,n2)).*exp(j*ep(:,n2)).*exp(j*2*pi*d/lambda*((0:M-1)')*(sin(phi_el)-sin(phi0)))); %(3.2.20) 
            for k1=1:K
                for k2=1:K
                    Rc_el(n1+N*(k1-1),n2+N*(k2-1))=exp(-Br^2*(k2-k1)^2/8)*(1+ca(n1,n2))*exp(j*cp(n1,n2))*...
                        G(n1)*conj(G(n2))*sum((abs(F).^2).*exp(j*(n1-n2)*omega_s+j*(k1-k2)*omega_d))*theta_step/(Rel^4);
                end
            end
        end
    end
    Rc=Rc+Rc_el;
end
Rc=Rc/max(max(abs(Rc)));R=Rc*(10^(CNR/10))+eye(N*K);

D=eig(R);D=sort(abs(D),'descend');D_db=10*log10(abs(D)/max(abs(D))); 
figure;plot(D_db,'*-'),grid on;
title('杂波特征谱');xlabel('特征数目');ylabel('特征值/dB');

psi_num=40;
fd_num=40;
psi_n=[-1:2/psi_num:1];
fd_n=[-1:2/fd_num:1];
R_inv=inv(R); 
for p=1:length(psi_n)
    Ss=exp(j*(0:N-1).'*pi*psi_n(p));
    for q=1:length(fd_n)
       St=exp(j*(0:K-1).'*pi*fd_n(q)); 
       S=kron(St,Ss);
       P(p,q)=1/abs(S'*R_inv*S);
    end
end
P_db=10*log10(P/max(max(P)));                     
% P_db(find(P_db<-60))=-60;                                       
figure;mesh(psi_n,fd_n,P_db); 
title('杂波功率谱');xlabel('2fd/fr');ylabel('cos(\psi)');zlabel('P/dB'); 
figure;contour(P_db,20);         
    
    
    