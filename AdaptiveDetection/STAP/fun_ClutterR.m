function R=fun_ClutterR(V,H,N,M,K,lambda,d,fr,phi0,theta0,Inn,Imm,Ikk,CNR,senser_error,channel_error,ea,ep,crab_theta,Br)
% fun_ClutterR is to generate covariance of the clutter.

V=140;                                            % �ػ��ٶ�
H=8000;                                           % �ػ��߶�
N=16;                                             % ÿ����N����Ԫ
M=8;                                              % ÿ����M����Ԫ
K=16;                                             % ʱ�������(ʱ��������)
lambda=0.23;                                      % ����
d=1/2*lambda;                                     % ��Ԫ����벨��֮��
c=3e+8;                                           % ����
fr=2434.8;                                        % PRF
Rmax=400000;                                      % �״�������þ���
phi0=0;                                           % �������ĸߵͽ�
theta0=pi/2;                                      % �������ķ�λ��
Inn=30;                                           % �����Ȩϵ����˥��dB��
Imm=30;                                           % �����Ȩϵ����˥��dB��
Ikk=30;                                           % ʱ���Ȩϵ����˥��dB��
CNR=50;                                           % �����

senser_error=0.0025;                              % ��Ԫ���ٷֱ�
channel_error=0.003;                              % ͨ�����ٷֱ�
ea=0.0025*randn(M,N);                             % ��Ԫ�������
ep=0.0005*randn(M,N);                             % ��Ԫ��λ���
crab_theta=0.1*pi/180;                            % �ɻ�ƫ����
Br=0.005;                                         % �Ӳ������
psi=0:pi/180:pi;                                  % ׶��ȡֵ��Χ
phi=0:pi/180:pi/2;                                % �ߵͽ�ȡֵ��Χ
theta=0:pi/180:pi;                                % ��λ��ȡֵ��Χ
theta_step=theta(2)-theta(1);
caa=(channel_error*randn(N,1));ca=caa*caa.';                 %ͨ����������(�м�������)
ca=ca-diag(diag(ca));
cpp=(channel_error*pi*randn(N,1));cp=cpp*cpp.';              %ͨ������λ���(�м���λ���)
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
    Rel=1.5*H+el*Ru;   % ���뻷����ʼλ�õ�ȡֵ��ΧΪ��[H,H+Ru]���˴���Ϊ��1.5H��
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
title('�Ӳ�������');xlabel('������Ŀ');ylabel('����ֵ/dB');

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
title('�Ӳ�������');xlabel('2fd/fr');ylabel('cos(\psi)');zlabel('P/dB'); 
figure;contour(P_db,20);         
    
    
    